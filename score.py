"""Score GRCh38 variants across AG2 reference conditions and scorers."""

import argparse
import csv
from datetime import datetime, timezone
from importlib.metadata import version
import io
from itertools import product
import json
from pathlib import Path
import sys
import tomllib

import grpc
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

from pilot import MODEL_VERSION, ORGANISM, api_key_from_file_or_environment
from reference import IndexedReference
from splicing import SPLICE_OUTPUTS, observed_summary


WINDOW_BP = dna_client.SEQUENCE_LENGTH_1MB
INPUT_COLUMNS = ("input_id", "chromosome", "position", "ref", "alt")
SCORE_IDS = ("input_id", "requested_condition", "requested_scorer")
REQUEST_COLUMNS = (*SCORE_IDS, "scored_interval", "status", "score_rows", "error")
SUMMARY_COLUMNS = (
    "input_id", "requested_condition", *SPLICE_OUTPUTS,
    "merged_splicing", "missing_observed_components",
)


def read_variants(
    text: str, reference: IndexedReference
) -> tuple[dict[str, genome.Variant], list[tuple[int, str, str]]]:
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
    headers = reader.fieldnames or []
    if set(INPUT_COLUMNS).difference(headers) or len(set(headers)) != len(headers):
        raise ValueError(f"Input TSV requires unique headers including {', '.join(INPUT_COLUMNS)}.")

    variants = {}
    rejected = []
    seen_ids = set()
    for row in reader:
        input_id = (row["input_id"] or "").strip()
        try:
            if None in row or any(row[name] is None for name in INPUT_COLUMNS):
                raise ValueError("Row does not match the TSV header.")
            if not input_id:
                raise ValueError("input_id is empty.")
            if input_id in seen_ids:
                raise ValueError("Duplicate input_id.")
            seen_ids.add(input_id)
            chromosome = row["chromosome"].strip()
            position = int(row["position"])
            ref = row["ref"].strip().upper()
            alt = row["alt"].strip().upper()
            if not chromosome or position < 1:
                raise ValueError("Expected a chromosome and a positive 1-based position.")
            if not ref or not alt or ref == alt or set(ref + alt).difference("ACGT"):
                raise ValueError("Expected distinct, nonempty REF and ALT sequences using A/C/G/T.")
            if len(ref) > WINDOW_BP:
                raise ValueError("REF does not fit in the scoring window.")
            variant = genome.Variant(
                chromosome=chromosome, position=position, reference_bases=ref, alternate_bases=alt
            )
            reference.verify_variant(variant)
        except ValueError as exc:
            rejected.append((reader.line_num, input_id, str(exc)))
        else:
            variants[input_id] = variant
    if not variants and not rejected:
        raise ValueError("Input TSV contains no variant rows.")
    return variants, rejected


def read_panel(
    path: Path, supported_scorers: list[str]
) -> tuple[str, int | str, dict[str, list[str]]]:
    with path.open("rb") as handle:
        panel = tomllib.load(handle)

    name = panel.get("name")
    panel_version = panel.get("version")
    conditions = panel.get("conditions")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Panel requires a nonempty name.")
    if type(panel_version) is not int and not (
        isinstance(panel_version, str) and panel_version.strip()
    ):
        raise ValueError("Panel version must be an integer or a nonempty string.")
    if not isinstance(conditions, dict) or not conditions:
        raise ValueError("Panel requires a nonempty [conditions] table.")
    for condition, names in conditions.items():
        if not condition or condition != condition.strip():
            raise ValueError("Panel condition IDs must be nonempty and have no surrounding spaces.")
        if not isinstance(names, list) or not names or not all(isinstance(n, str) for n in names):
            raise ValueError(f"Panel condition {condition} requires a nonempty list of scorer names.")
        if len(set(names)) != len(names):
            raise ValueError(f"Panel condition {condition} contains duplicate scorers.")
        unknown = set(names).difference(supported_scorers)
        if unknown:
            raise ValueError(f"Unsupported scorer(s) for {condition}: {', '.join(sorted(unknown))}.")
    return name, panel_version, conditions


def run_batch(
    args: argparse.Namespace,
    condition_scorers: dict[str, list[str]],
    panel: dict | None,
) -> bool:
    input_text = args.input.read_text(encoding="utf-8-sig")
    # Explicit SV alleles can exceed the CSV parser's default field limit.
    csv.field_size_limit(max(csv.field_size_limit(), len(input_text)))
    reference = IndexedReference(args.reference)
    variants, rejected = read_variants(input_text, reference)
    planned_calls = len(variants) * len(condition_scorers)
    if args.dry_run:
        print(f"REF-verified rows: {len(variants)}")
        print(f"Rejected rows: {len(rejected)}")
        print(f"Conditions: {len(condition_scorers)}")
        print(f"Planned score calls: {planned_calls}")
        return not rejected

    recommended = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    output = args.output_dir
    if output.exists():
        raise ValueError(f"Choose a new output directory: {output}")

    if variants:
        model = dna_client.create(
            api_key_from_file_or_environment(args.api_key_file), model_version=MODEL_VERSION, timeout=30
        )
        catalog = model.supported_conditions(ORGANISM).set_index("condition")
        unknown = set(condition_scorers).difference(catalog.index)
        if unknown:
            raise ValueError(f"Unknown AG2 reference conditions: {', '.join(sorted(unknown))}")
        conditions = catalog.loc[list(condition_scorers)]

    output.mkdir(parents=True)
    (output / "variants.tsv").write_text(input_text, encoding="utf-8")
    if variants:
        conditions.to_csv(output / "conditions.tsv", sep="\t", encoding="utf-8")
    if rejected:
        pd.DataFrame(rejected, columns=["line", "input_id", "error"]).to_csv(
            output / "rejected.tsv", sep="\t", index=False, encoding="utf-8"
        )
        print(f"Rejected {len(rejected)} input rows. See rejected.tsv.", file=sys.stderr)

    requests = []
    tables = []
    splicing_conditions = {
        condition for condition, names in condition_scorers.items()
        if set(SPLICE_OUTPUTS).issubset(names)
    }
    summaries = []
    api_error = None
    attempted_calls = 0
    started_at = datetime.now(timezone.utc)
    for (input_id, variant), (condition, names) in product(variants.items(), condition_scorers.items()):
        interval = variant.reference_interval.resize(WINDOW_BP)
        scorers = [recommended[name] for name in names]
        attempted_calls += 1
        try:
            response = model.score_variant(
                interval=interval, variant=variant, variant_scorers=scorers,
                organism=ORGANISM, condition=condition,
            )
        except grpc.RpcError as exc:
            api_error = exc.code().name
            requests.extend(
                (input_id, condition, name, str(interval), "api_error", 0, api_error)
                for name in names
            )
            print(f"{input_id}, {condition}: API returned {api_error}. Stopped.", file=sys.stderr)
            break

        pair_tables = []
        for name, result in zip(names, response, strict=True):
            scores = variant_scorers.tidy_anndata(result, match_gene_strand=True)
            status = "no_rows" if scores.empty else "scored"
            requests.append((input_id, condition, name, str(interval), status, len(scores), ""))
            if not scores.empty:
                scores.insert(0, "requested_scorer", name)
                scores.insert(0, "requested_condition", condition)
                scores.insert(0, "input_id", input_id)
                pair_tables.append(scores)
        pair_scores = pd.concat(pair_tables, ignore_index=True) if pair_tables else pd.DataFrame()
        if not pair_scores.empty:
            tables.append(pair_scores)
        if condition in splicing_conditions:
            summaries.append({
                "input_id": input_id, "requested_condition": condition,
                **observed_summary(pair_scores),
            })
        print(f"{input_id}, {condition}: scored {len(scorers)} requested scorers", flush=True)

    scores = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(
        columns=[*SCORE_IDS, "raw_score", "imputation_mode"]
    )
    scores.to_csv(output / "scores.tsv", sep="\t", index=False, encoding="utf-8")
    pd.DataFrame(requests, columns=REQUEST_COLUMNS).to_csv(
        output / "requests.tsv", sep="\t", index=False, encoding="utf-8"
    )
    if splicing_conditions:
        pd.DataFrame(summaries, columns=SUMMARY_COLUMNS).to_csv(
            output / "splicing_summary.tsv", sep="\t", index=False, encoding="utf-8"
        )
    metadata = {
        "started_at_utc": started_at.isoformat(),
        "assembly": "GRCh38",
        "organism": ORGANISM.name,
        "model_version": MODEL_VERSION.name,
        "sdk_version": version("alphagenome"),
        "window_bp": WINDOW_BP,
        "conditions": list(condition_scorers),
        "scorers": list(dict.fromkeys(name for names in condition_scorers.values() for name in names)),
        "condition_scorers": condition_scorers,
        "splicing_summary_tracks": "observed" if splicing_conditions else None,
        "verified_variants": len(variants),
        "rejected_rows": len(rejected),
        "planned_score_calls": planned_calls,
        "attempted_score_calls": attempted_calls,
        "api_error": api_error,
        "match_gene_strand": True,
    }
    if panel is not None:
        metadata["panel"] = panel
    manifest = output / "request.json.tmp"
    manifest.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    manifest.rename(output / "request.json")
    return not (rejected or api_error)


def main() -> int:
    supported = variant_scorers.get_recommended_scorers(ORGANISM.to_proto())
    choices = [
        name for name, scorer in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items()
        if scorer in supported
    ]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="UTF-8 variant TSV")
    parser.add_argument("--reference", type=Path, required=True, help="Indexed GRCh38 FASTA")
    parser.add_argument("--condition", action="append", help="AG2 reference condition; repeatable")
    parser.add_argument("--scorer", action="append", choices=choices, help="AG2 scorer; repeatable")
    parser.add_argument("--panel", type=Path, help="TOML mapping of reference conditions to scorers")
    parser.add_argument("--dry-run", action="store_true", help="Verify REF and count planned score calls offline")
    parser.add_argument("--output-dir", type=Path, help="New directory for results; required for scoring")
    parser.add_argument("--api-key-file", type=Path, help="Private text file containing the API key")
    args = parser.parse_args()
    if args.panel is None:
        if not args.condition or not args.scorer:
            parser.error("Specify --panel or both --condition and --scorer.")
        for name in ("condition", "scorer"):
            values = getattr(args, name)
            if len(set(values)) != len(values):
                parser.error(f"Duplicate --{name} values.")
    elif args.condition or args.scorer:
        parser.error("--panel cannot be used with --condition or --scorer.")
    if not args.dry_run and args.output_dir is None:
        parser.error("--output-dir is required for scoring.")

    try:
        panel = None
        if args.panel is None:
            condition_scorers = dict.fromkeys(args.condition, args.scorer)
        else:
            name, panel_version, condition_scorers = read_panel(args.panel, choices)
            panel = {"name": name, "version": panel_version, "conditions": condition_scorers}
        complete = run_batch(args, condition_scorers, panel)
    except grpc.FutureTimeoutError:
        print("Timed out connecting to the AlphaGenome API.", file=sys.stderr)
        return 1
    except grpc.RpcError as exc:
        print(f"AlphaGenome API returned {exc.code().name}.", file=sys.stderr)
        return 1
    except (OSError, ValueError, csv.Error) as exc:
        print(exc, file=sys.stderr)
        return 1
    if not args.dry_run:
        print(f"Results saved to {args.output_dir}")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
