"""Score one public GRCh38 SNV in two AlphaGenome 2 reference conditions."""

import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import sys

import grpc
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

from reference import verify_snv_reference


VARIANT = genome.Variant(
    chromosome="chr2", position=219243701, reference_bases="C", alternate_bases="T"
)
TARGET_GENE = "ENSG00000163521"  # GLB1L
CONDITIONS = ("UBERON:0000178_CL:0001054", "CL:0000115")
SCORER = "RNA_SEQ"
MODEL_VERSION = dna_client.ModelVersion.ALL_FOLDS
ORGANISM = dna_client.Organism.HOMO_SAPIENS


def api_key_from_file_or_environment(path: Path | None) -> str:
    if path is not None:
        key = path.read_text(encoding="utf-8-sig").strip()
    else:
        key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        raise ValueError("Set ALPHAGENOME_API_KEY or pass a private --api-key-file.")
    return key


def run_pilot(reference: Path, api_key_file: Path | None, output: Path) -> None:
    if output.exists():
        raise ValueError(f"Choose a new output directory: {output}")
    verify_snv_reference(reference, VARIANT)

    interval = VARIANT.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    model = dna_client.create(
        api_key_from_file_or_environment(api_key_file), model_version=MODEL_VERSION, timeout=30
    )
    catalog = model.supported_conditions(ORGANISM)
    selected = catalog.set_index("condition").reindex(list(CONDITIONS))[
        ["biosample_name", "rna_seq"]
    ]
    if not selected["rna_seq"].eq("observed").fillna(False).all():
        raise ValueError("Both pilot conditions must have an observed RNA_SEQ track.")

    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[SCORER]
    sdk_version = version("alphagenome")
    tables = []
    for condition in CONDITIONS:
        response = model.score_variant(
            interval=interval,
            variant=VARIANT,
            variant_scorers=[scorer],
            organism=ORGANISM,
            condition=condition,
        )
        scores = variant_scorers.tidy_scores(response, match_gene_strand=True)
        if scores is None or scores.empty:
            raise ValueError(f"Empty RNA_SEQ response for {condition}.")

        observed_target = scores.loc[
            scores["gene_id"].eq(TARGET_GENE) & scores["imputation_mode"].eq("observed"),
            "raw_score",
        ]
        if not observed_target.dropna().map(math.isfinite).any():
            raise ValueError(f"No finite GLB1L score on an observed track for {condition}.")

        scores.insert(0, "requested_condition", condition)
        scores.insert(1, "sdk_version", sdk_version)
        tables.append(scores)
        print(f"{condition}: {len(scores)} score rows", flush=True)

    metadata = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "assembly": "GRCh38",
        "variant": str(VARIANT),
        "reference_allele_verified": VARIANT.reference_bases,
        "interval_0_based_half_open": [interval.chromosome, interval.start, interval.end],
        "target_gene_id": TARGET_GENE,
        "conditions": list(CONDITIONS),
        "scorer": SCORER,
        "model_version": MODEL_VERSION.name,
        "organism": ORGANISM.name,
        "sdk_version": sdk_version,
        "python_version": platform.python_version(),
        "match_gene_strand": True,
    }

    output.mkdir(parents=True)
    selected.to_csv(output / "conditions.tsv", sep="\t", encoding="utf-8")
    pd.concat(tables, ignore_index=True).to_csv(
        output / "scores.tsv", sep="\t", index=False, encoding="utf-8"
    )
    # Write the manifest last so its presence marks a complete run.
    manifest = output / "request.json.tmp"
    manifest.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    manifest.rename(output / "request.json")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="Indexed GRCh38 FASTA")
    parser.add_argument("--api-key-file", type=Path, help="Private text file containing the API key")
    parser.add_argument("--output-dir", required=True, type=Path, help="New directory for results")
    args = parser.parse_args()
    try:
        run_pilot(args.reference, args.api_key_file, args.output_dir)
    except grpc.FutureTimeoutError:
        print("Timed out connecting to the AlphaGenome API.", file=sys.stderr)
        return 1
    except grpc.RpcError as exc:
        print(f"AlphaGenome API returned {exc.code().name}.", file=sys.stderr)
        return 1
    except (OSError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"Results saved to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
