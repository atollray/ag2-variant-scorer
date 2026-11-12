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

import grpc
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers

from pilot import MODEL_VERSION, ORGANISM, api_key_from_file_or_environment
from reference import IndexedReference


WINDOW_BP = dna_client.SEQUENCE_LENGTH_1MB
INPUT_COLUMNS = ("input_id", "chromosome", "position", "ref", "alt")
SCORE_IDS = ("input_id", "requested_condition", "requested_scorer")
REQUEST_COLUMNS = (*SCORE_IDS, "scored_interval", "status", "score_rows", "error")


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


def run_batch(args: argparse.Namespace) -> bool:
    input_text = args.input.read_text(encoding="utf-8-sig")
    # Explicit SV alleles can exceed the CSV parser's default field limit.
    csv.field_size_limit(max(csv.field_size_limit(), len(input_text)))
    reference = IndexedReference(args.reference)
    variants, rejected = read_variants(input_text, reference)
    recommended = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    scorers = [recommended[name] for name in args.scorer]
    output = args.output_dir
    if output.exists():
        raise ValueError(f"Choose a new output directory: {output}")

    if variants:
        model = dna_client.create(
            api_key_from_file_or_environment(args.api_key_file), model_version=MODEL_VERSION, timeout=30
        )
        catalog = model.supported_conditions(ORGANISM).set_index("condition")
        unknown = set(args.condition).difference(catalog.index)
        if unknown:
            raise ValueError(f"Unknown AG2 reference conditions: {', '.join(sorted(unknown))}")
        conditions = catalog.loc[args.condition]

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
    api_error = None
    started_at = datetime.now(timezone.utc)
    for (input_id, variant), condition in product(variants.items(), args.condition):
        interval = variant.reference_interval.resize(WINDOW_BP)
        try:
            response = model.score_variant(
                interval=interval, variant=variant, variant_scorers=scorers,
                organism=ORGANISM, condition=condition,
            )
        except grpc.RpcError as exc:
            api_error = exc.code().name
            requests.extend(
                (input_id, condition, name, str(interval), "api_error", 0, api_error)
                for name in args.scorer
            )
            print(f"{input_id}, {condition}: API returned {api_error}. Stopped.", file=sys.stderr)
            break

        for name, result in zip(args.scorer, response, strict=True):
            scores = variant_scorers.tidy_anndata(result, match_gene_strand=True)
            status = "no_rows" if scores.empty else "scored"
            requests.append((input_id, condition, name, str(interval), status, len(scores), ""))
            if not scores.empty:
                scores.insert(0, "requested_scorer", name)
                scores.insert(0, "requested_condition", condition)
                scores.insert(0, "input_id", input_id)
                tables.append(scores)
        print(f"{input_id}, {condition}: scored {len(scorers)} requested scorers", flush=True)

    scores = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(
        columns=[*SCORE_IDS, "raw_score", "imputation_mode"]
    )
    scores.to_csv(output / "scores.tsv", sep="\t", index=False, encoding="utf-8")
    pd.DataFrame(requests, columns=REQUEST_COLUMNS).to_csv(
        output / "requests.tsv", sep="\t", index=False, encoding="utf-8"
    )
    metadata = {
        "started_at_utc": started_at.isoformat(),
        "assembly": "GRCh38",
        "organism": ORGANISM.name,
        "model_version": MODEL_VERSION.name,
        "sdk_version": version("alphagenome"),
        "window_bp": WINDOW_BP,
        "conditions": args.condition,
        "scorers": args.scorer,
        "verified_variants": len(variants),
        "rejected_rows": len(rejected),
        "planned_score_calls": len(variants) * len(args.condition),
        "attempted_score_calls": len(requests) // len(scorers),
        "api_error": api_error,
        "match_gene_strand": True,
    }
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
    parser.add_argument("--condition", action="append", required=True, help="AG2 reference condition; repeatable")
    parser.add_argument("--scorer", action="append", choices=choices, required=True, help="AG2 scorer; repeatable")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory for results")
    parser.add_argument("--api-key-file", type=Path, help="Private text file containing the API key")
    args = parser.parse_args()
    for name in ("condition", "scorer"):
        values = getattr(args, name)
        if len(set(values)) != len(values):
            parser.error(f"Duplicate --{name} values.")
    try:
        complete = run_batch(args)
    except grpc.FutureTimeoutError:
        print("Timed out connecting to the AlphaGenome API.", file=sys.stderr)
        return 1
    except grpc.RpcError as exc:
        print(f"AlphaGenome API returned {exc.code().name}.", file=sys.stderr)
        return 1
    except (OSError, ValueError, csv.Error) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"Results saved to {args.output_dir}")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
