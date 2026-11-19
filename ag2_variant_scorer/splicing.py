"""Run the AG2 splicing pilot and save scores and per-allele predictions."""

from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

from .reference import IndexedReference
from .sv_examples import AFFECTED_REGION, prepare_sv_examples


SPLICE_OUTPUTS = ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
CONDITION = "UBERON:0001159"  # Sigmoid colon in the official APOL4 example.
CONTROL = genome.Variant(
    chromosome="chr22", position=36201698, reference_bases="A", alternate_bases="C"
)


def splice_summary(maxima) -> dict:
    components = {name: maxima.get(name) for name in SPLICE_OUTPUTS}
    missing = [name for name, value in components.items() if value is None]
    merged = None
    if not missing:
        sites, usage, junctions = (components[name] for name in SPLICE_OUTPUTS)
        merged = sites + usage + junctions / 5
    return {
        **components,
        "merged_splicing": merged,
        "missing_observed_components": ",".join(missing),
    }


def observed_summary(scores: pd.DataFrame) -> dict:
    maxima = {}
    if not scores.empty:
        observed = scores.loc[scores["imputation_mode"].eq("observed")]
        observed = observed.loc[np.isfinite(observed["raw_score"])]
        maxima = observed["raw_score"].abs().groupby(observed["output_type"]).max().dropna()
    return splice_summary(maxima)


def junction_row(interval: genome.Interval) -> tuple:
    """Read the SDK's junction boundaries from an Interval and its strand."""
    donor = acceptor = None
    if interval.strand == "+":
        donor, acceptor = interval.start, interval.end
    elif interval.strand == "-":
        donor, acceptor = interval.end, interval.start
    return interval.chromosome, interval.start, interval.end, interval.strand, donor, acceptor


def save_allele_predictions(path: Path, prediction) -> dict:
    """Keep values, track labels and SDK-reported coordinates in matching order."""
    arrays = {}
    geometry = {}
    for name in SPLICE_OUTPUTS:
        key = name.lower()
        data = getattr(prediction, key)
        if data is None:
            raise ValueError(f"Missing prediction for {path.name} {name}.")
        values = np.asarray(data.values)
        if values.dtype.name == "bfloat16":
            values = values.astype(np.float32)
        arrays[key] = values
        data.metadata.to_csv(path / f"{key}_tracks.tsv", sep="\t", index=False, encoding="utf-8")

        interval = data.interval
        geometry[key] = {
            "interval_0_based_half_open": [interval.chromosome, interval.start, interval.end]
            if interval is not None else None,
            "interval_strand": interval.strand if interval is not None else None,
            "shape": list(values.shape),
            "dtype": str(values.dtype),
        }
        if name == "SPLICE_JUNCTIONS":
            coordinates = pd.DataFrame.from_records(
                (junction_row(junction) for junction in data.junctions),
                columns=["chromosome", "start0", "end0", "strand", "donor_boundary0", "acceptor_boundary0"],
            )
            coordinates.to_csv(path / "junctions.tsv", sep="\t", index=False, encoding="utf-8")
        else:
            geometry[key]["resolution_bp"] = data.resolution

    np.savez_compressed(path / "values.npz", **arrays)
    return geometry


def run_splicing_pilot(
    fasta: Path, api_key: str, output: Path, *, include_sv: bool = False
) -> None:
    if output.exists():
        raise ValueError(f"Choose a new output directory: {output}")
    reference = IndexedReference(fasta)
    reference.verify_variant(CONTROL)
    cases = {"APOL4_SNV": CONTROL}
    if include_sv:
        cases.update(prepare_sv_examples(reference))

    organism = dna_client.Organism.HOMO_SAPIENS
    model_version = dna_client.ModelVersion.ALL_FOLDS
    model = dna_client.create(api_key, model_version=model_version, timeout=30)
    catalog = model.supported_conditions(organism).set_index("condition")
    condition_metadata = catalog.reindex([CONDITION])[
        ["biosample_name", *[name.lower() for name in SPLICE_OUTPUTS]]
    ]
    if not condition_metadata.iloc[0, 1:].eq("observed").fillna(False).all():
        raise ValueError("The APOL4 pilot requires all three observed splicing outputs in sigmoid colon.")

    output.mkdir(parents=True)
    condition_metadata.to_csv(output / "conditions.tsv", sep="\t", encoding="utf-8")
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in SPLICE_OUTPUTS]
    requested_outputs = [dna_client.OutputType[name] for name in SPLICE_OUTPUTS]
    summaries = []
    request_cases = {}
    sdk_version = version("alphagenome")
    for case_id, variant in cases.items():
        center = AFFECTED_REGION if case_id != "APOL4_SNV" else variant.reference_interval
        interval = center.resize(dna_client.SEQUENCE_LENGTH_1MB)
        case_dir = output / case_id
        case_dir.mkdir()
        scores = variant_scorers.tidy_scores(
            model.score_variant(
                interval=interval, variant=variant, variant_scorers=scorers,
                condition=CONDITION, organism=organism,
            ),
            match_gene_strand=True,
        )
        if scores is None or scores.empty:
            raise ValueError(f"Empty splicing score response for {case_id}.")
        scores.insert(0, "variant_case", case_id)
        scores.insert(1, "requested_condition", CONDITION)
        scores.insert(2, "sdk_version", sdk_version)
        scores.to_csv(case_dir / "scores.tsv", sep="\t", index=False, encoding="utf-8")
        summaries.append({
            "variant_case": case_id, "requested_condition": CONDITION, **observed_summary(scores)
        })

        prediction = model.predict_variant(
            interval=interval, variant=variant, requested_outputs=requested_outputs,
            condition=CONDITION, organism=organism,
        )
        allele_geometry = {}
        for allele, result in (("REF", prediction.reference), ("ALT", prediction.alternate)):
            path = case_dir / allele
            path.mkdir()
            allele_geometry[allele] = save_allele_predictions(path, result)

        request_cases[case_id] = {
            "chromosome": variant.chromosome,
            "position_1_based": variant.position,
            "reference_bases": variant.reference_bases,
            "alternate_bases": variant.alternate_bases,
            "interval_0_based_half_open": [interval.chromosome, interval.start, interval.end],
            "sdk_prediction_geometry_by_allele": allele_geometry,
        }
        print(f"{case_id}: saved {len(scores)} score rows and REF/ALT predictions", flush=True)

    pd.DataFrame(summaries).to_csv(output / "summary.tsv", sep="\t", index=False, encoding="utf-8")
    metadata = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "assembly": "GRCh38",
        "condition": CONDITION,
        "organism": organism.name,
        "model_version": model_version.name,
        "sdk_version": sdk_version,
        "scorers_and_outputs": list(SPLICE_OUTPUTS),
        "match_gene_strand": True,
        "summary_tracks": "observed",
        "prediction_coordinates": "SDK-reported per allele; no cross-allele remapping performed",
        "cases": request_cases,
    }
    if include_sv:
        metadata["sv_affected_region_1_based_inclusive"] = [
            AFFECTED_REGION.chromosome, AFFECTED_REGION.start + 1, AFFECTED_REGION.end
        ]
    manifest = output / "request.json.tmp"
    manifest.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    manifest.rename(output / "request.json")
