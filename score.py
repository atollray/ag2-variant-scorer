"""Score GRCh38 variants across AG2 reference conditions and scorers."""

import argparse
from contextlib import closing, ExitStack, nullcontext
import csv
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import io
import json
from pathlib import Path
import sqlite3
import sys
import tomllib

import grpc
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, variant_scorers
from filelock import Timeout

from pilot import MODEL_VERSION, ORGANISM, api_key_from_file_or_environment
from reference import IndexedReference
from splicing import SPLICE_OUTPUTS, observed_summary
from storage import CachedScore, Outcome, RunStore, SharedCache, run_lock


WINDOW_BP = dna_client.SEQUENCE_LENGTH_1MB
INPUT_COLUMNS = ("input_id", "chromosome", "position", "ref", "alt")
SCORE_IDS = ("input_id", "requested_condition", "requested_scorer")
REQUEST_COLUMNS = (*SCORE_IDS, "scored_interval", "status", "score_rows", "error", "source", "retrieved_at_utc")
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


def scorer_fingerprints(condition_scorers: dict[str, list[str]]) -> dict[str, str]:
    recommended = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    names = dict.fromkeys(name for group in condition_scorers.values() for name in group)
    return {
        name: recommended[name].to_proto().SerializeToString(deterministic=True).hex()
        for name in names
    }


def scoring_settings() -> dict:
    return {
        "assembly": "GRCh38",
        "organism": ORGANISM.name,
        "model_version": MODEL_VERSION.name,
        "sdk_version": version("alphagenome"),
        "window_bp": WINDOW_BP,
        "match_gene_strand": True,
    }


def work_items(plan: dict):
    metadata = plan["metadata"]
    for input_id, chromosome, position, ref, alt in plan["variants"]:
        variant = genome.Variant(
            chromosome=chromosome, position=position, reference_bases=ref, alternate_bases=alt
        )
        interval = variant.reference_interval.resize(metadata["window_bp"])
        for condition, names in metadata["condition_scorers"].items():
            yield input_id, variant, interval, condition, names


def request_key(variant, interval, condition: str, settings: dict) -> str:
    identity = {
        **settings,
        "format": "ag2-tidy-tsv-v1",
        "variant": [
            variant.chromosome, variant.position, variant.reference_bases, variant.alternate_bases,
        ],
        "interval": [interval.chromosome, interval.start, interval.end, interval.strand],
        "condition": condition,
    }
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def scorer_key(plan: dict, name: str) -> str:
    return f"{name}:{plan['scorer_fingerprints'][name]}"


def preview(plan: dict, cache_dir: Path | None) -> None:
    metadata = plan["metadata"]
    print(f"REF-verified rows: {metadata['verified_variants']}")
    print(f"Rejected rows: {metadata['rejected_rows']}")
    print(f"Conditions: {len(metadata['conditions'])}")
    print(f"Planned score calls: {metadata['planned_score_calls']}")
    if cache_dir is None:
        return
    cached = uncached = 0
    settings = scoring_settings()
    with closing(SharedCache(cache_dir, readonly=True)) as cache:
        for _, variant, interval, condition, names in work_items(plan):
            key = request_key(variant, interval, condition, settings)
            hits = sum(cache.contains(key, scorer_key(plan, name)) for name in names)
            cached += hits
            uncached += hits < len(names)
    print(f"Cached scorer results: {cached}")
    print(f"Uncached score calls (upper bound): {uncached}")


def condition_metadata(model, conditions: list[str]) -> dict[str, str]:
    catalog = model.supported_conditions(ORGANISM).set_index("condition")
    unknown = set(conditions).difference(catalog.index)
    if unknown:
        raise ValueError(f"Unknown AG2 reference conditions: {', '.join(sorted(unknown))}")
    return {
        condition: catalog.loc[[condition]].to_csv(sep="\t", lineterminator="\n")
        for condition in conditions
    }


def save_cache(cache: SharedCache, group_key: str, plan: dict, condition_tsv: str, outcomes: dict) -> None:
    entries = {
        scorer_key(plan, name): CachedScore(
            result.payload, result.row_count, condition_tsv, result.retrieved_at_utc
        )
        for name, result in outcomes.items() if result.status == "scored"
    }
    if entries:
        cache.save(group_key, entries)


def acquire_scores(
    state: RunStore, plan: dict, cache: SharedCache | None, api_key_file: Path | None, *, restore_cache: bool
):
    model = None
    settings = scoring_settings()
    condition_rows = state.get_conditions()
    conditions = plan["metadata"]["conditions"]
    recommended = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    for input_id, variant, interval, condition, names in work_items(plan):
        previous = state.get_group(input_id, condition)
        pending = [name for name in names if name not in previous or previous[name].status == "api_error"]
        if not pending and (cache is None or not restore_cache):
            continue

        group_key = request_key(variant, interval, condition, settings) if cache is not None else ""
        backfill = {}
        if cache is not None and restore_cache:
            backfill = {
                name: result for name, result in previous.items()
                if result.status == "scored" and result.source == "api"
                and not cache.contains(group_key, scorer_key(plan, name))
            }
        if not pending and not backfill:
            continue

        with cache.lock(group_key) if cache else nullcontext():
            hits = {}
            if cache is not None:
                if backfill:
                    # The run can have committed its response just before cache publication stopped.
                    save_cache(cache, group_key, plan, condition_rows[condition], backfill)
                for name in pending:
                    hit = cache.lookup(group_key, scorer_key(plan, name))
                    if hit is not None:
                        hits[name] = Outcome(
                            "scored", hit.payload, hit.row_count, source="cache",
                            retrieved_at_utc=hit.retrieved_at_utc,
                        )
                        if condition not in condition_rows:
                            condition_rows[condition] = hit.condition_tsv
                            state.set_conditions({condition: hit.condition_tsv})
                if hits:
                    state.record(input_id, condition, hits)
                    pending = [name for name in pending if name not in hits]
            if not pending:
                continue

            if model is None:
                model = dna_client.create(
                    api_key_from_file_or_environment(api_key_file), model_version=MODEL_VERSION, timeout=30
                )
                if set(conditions).difference(condition_rows):
                    condition_rows = condition_metadata(model, conditions)
                    state.set_conditions(condition_rows)

            state.start_attempt()
            try:
                response = model.score_variant(
                    interval=interval, variant=variant,
                    variant_scorers=[recommended[name] for name in pending],
                    organism=ORGANISM, condition=condition,
                )
            except grpc.RpcError as exc:
                code = exc.code().name
                state.record(input_id, condition, {
                    name: Outcome("api_error", error=code, source="api") for name in pending
                })
                print(f"{input_id}, {condition}: API returned {code}. Stopped.", file=sys.stderr)
                return code

            outcomes = {}
            retrieved_at = datetime.now(timezone.utc).isoformat()
            for name, result in zip(pending, response, strict=True):
                scores = variant_scorers.tidy_anndata(result, match_gene_strand=True)
                if scores.empty:
                    outcomes[name] = Outcome("no_rows", source="api", retrieved_at_utc=retrieved_at)
                else:
                    payload = scores.to_csv(sep="\t", index=False, lineterminator="\n")
                    outcomes[name] = Outcome(
                        "scored", payload, len(scores), source="api", retrieved_at_utc=retrieved_at
                    )
            state.record(input_id, condition, outcomes)
            if cache is not None:
                save_cache(cache, group_key, plan, condition_rows[condition], outcomes)
            print(f"{input_id}, {condition}: scored {len(pending)} requested scorers", flush=True)
    return None


def score_fields(state: RunStore, plan: dict) -> list[str]:
    fields = dict.fromkeys(SCORE_IDS)
    for input_id, _, _, condition, names in work_items(plan):
        group = state.get_group(input_id, condition)
        for name in names:
            outcome = group.get(name)
            if outcome is not None and outcome.payload:
                header = next(csv.reader(io.StringIO(outcome.payload), delimiter="\t"))
                fields.update(dict.fromkeys(header))
    if len(fields) == len(SCORE_IDS):
        fields.update(dict.fromkeys(("raw_score", "imputation_mode")))
    return list(fields)


def write_exports(output: Path, state: RunStore, plan: dict, error: str | None, interrupted: bool) -> None:
    metadata = dict(plan["metadata"])
    splicing_conditions = {
        condition for condition, names in metadata["condition_scorers"].items()
        if set(SPLICE_OUTPUTS).issubset(names)
    }
    headers = {"scores.tsv": score_fields(state, plan), "requests.tsv": REQUEST_COLUMNS}
    if splicing_conditions:
        headers["splicing_summary.tsv"] = SUMMARY_COLUMNS

    cached = 0
    with ExitStack() as files:
        writers = {}
        for filename, columns in headers.items():
            handle = files.enter_context((output / f"{filename}.tmp").open("w", encoding="utf-8", newline=""))
            writers[filename] = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writers[filename].writeheader()

        for input_id, _, interval, condition, names in work_items(plan):
            group = state.get_group(input_id, condition)
            splicing = []
            for name in names:
                outcome = group.get(name)
                if outcome is None:
                    continue
                identifiers = dict(zip(SCORE_IDS, (input_id, condition, name), strict=True))
                writers["requests.tsv"].writerow({
                    **identifiers, "scored_interval": str(interval), "status": outcome.status,
                    "score_rows": outcome.row_count, "error": outcome.error, "source": outcome.source,
                    "retrieved_at_utc": outcome.retrieved_at_utc,
                })
                cached += outcome.source == "cache"
                if outcome.payload is None:
                    continue
                csv.field_size_limit(max(csv.field_size_limit(), len(outcome.payload)))
                rows = csv.DictReader(io.StringIO(outcome.payload, newline=""), delimiter="\t")
                for row in rows:
                    writers["scores.tsv"].writerow({**row, **identifiers})
                    if condition in splicing_conditions and name in SPLICE_OUTPUTS:
                        splicing.append({
                            "output_type": row["output_type"],
                            "imputation_mode": row["imputation_mode"],
                            "raw_score": float(row["raw_score"] or "nan"),
                        })
            if condition in splicing_conditions and all(
                name in group and group[name].status != "api_error" for name in SPLICE_OUTPUTS
            ):
                writers["splicing_summary.tsv"].writerow({
                    "input_id": input_id, "requested_condition": condition,
                    **observed_summary(pd.DataFrame(splicing)),
                })

    tables = list(headers)
    if plan["rejected"]:
        pd.DataFrame(plan["rejected"], columns=["line", "input_id", "error"]).to_csv(
            output / "rejected.tsv.tmp", sep="\t", index=False, encoding="utf-8"
        )
        tables.append("rejected.tsv")
    condition_rows = state.get_conditions()
    if condition_rows:
        catalog = pd.concat([
            pd.read_csv(io.StringIO(condition_rows[condition]), sep="\t", dtype=str, keep_default_na=False)
            for condition in metadata["conditions"] if condition in condition_rows
        ], ignore_index=True)
        catalog.to_csv(output / "conditions.tsv.tmp", sep="\t", index=False, encoding="utf-8")
        tables.append("conditions.tsv")

    metadata.update(
        attempted_score_calls=state.attempted_calls, cached_scorer_results=cached,
        splicing_summary_tracks="observed" if splicing_conditions else None,
        api_error=error, interrupted=interrupted,
    )
    manifest = output / "request.json"
    pending_manifest = output / "request.json.tmp"
    pending_manifest.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    # The manifest marks a consistent set of exports. SQLite remains the resume source.
    manifest.unlink(missing_ok=True)
    for filename in tables:
        (output / f"{filename}.tmp").replace(output / filename)
    pending_manifest.replace(manifest)


def execute(
    output: Path, state: RunStore, plan: dict, cache_dir: Path | None, api_key_file: Path | None,
    *, restore_cache: bool = False,
) -> int:
    error = None
    interrupted = False
    try:
        open_cache = cache_dir is not None and bool(plan["variants"])
        with closing(SharedCache(cache_dir)) if open_cache else nullcontext() as cache:
            error = acquire_scores(state, plan, cache, api_key_file, restore_cache=restore_cache)
    except grpc.FutureTimeoutError:
        error = "CONNECT_TIMEOUT"
        print("Timed out connecting to the AlphaGenome API.", file=sys.stderr)
    except grpc.RpcError as exc:
        error = exc.code().name
        print(f"AlphaGenome API returned {error}.", file=sys.stderr)
    except KeyboardInterrupt:
        interrupted = True
        print("Interrupted. Resume with --resume and this output directory.", file=sys.stderr)

    write_exports(output, state, plan, error, interrupted)
    print(f"Results saved to {output}")
    if interrupted:
        return 130
    return 1 if plan["metadata"]["rejected_rows"] or error else 0


def run_batch(
    args: argparse.Namespace,
    condition_scorers: dict[str, list[str]],
    panel: dict | None,
) -> int:
    input_text = args.input.read_text(encoding="utf-8-sig")
    # Explicit SV alleles can exceed the CSV parser's default field limit.
    csv.field_size_limit(max(csv.field_size_limit(), len(input_text)))
    variants, rejected = read_variants(input_text, IndexedReference(args.reference))
    fingerprints = scorer_fingerprints(condition_scorers)
    metadata = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        **scoring_settings(),
        "conditions": list(condition_scorers),
        "scorers": list(fingerprints),
        "condition_scorers": condition_scorers,
        "verified_variants": len(variants),
        "rejected_rows": len(rejected),
        "planned_score_calls": len(variants) * len(condition_scorers),
    }
    if panel is not None:
        metadata["panel"] = panel
    plan = {
        "metadata": metadata,
        "variants": [
            [input_id, v.chromosome, v.position, v.reference_bases, v.alternate_bases]
            for input_id, v in variants.items()
        ],
        "input_sha256": hashlib.sha256(input_text.encode("utf-8")).hexdigest(),
        "scorer_fingerprints": fingerprints,
        "cache_dir": str(args.cache_dir.resolve()) if args.cache_dir else None,
        "rejected": rejected,
    }
    if args.dry_run:
        preview(plan, args.cache_dir)
        return int(bool(rejected))

    output = args.output_dir
    try:
        output.mkdir(parents=True)
    except FileExistsError as exc:
        raise ValueError(f"Choose a new output directory: {output}") from exc
    with run_lock(output):
        (output / "variants.tsv").write_text(input_text, encoding="utf-8")
        if rejected:
            print(f"Rejected {len(rejected)} input rows. See rejected.tsv.", file=sys.stderr)
        with closing(RunStore.create(output, plan)) as state:
            return execute(output, state, plan, args.cache_dir, args.api_key_file)


def resume(output: Path, cache_dir: Path | None, api_key_file: Path | None) -> int:
    if not output.is_dir():
        raise ValueError(f"Run directory does not exist: {output}")
    with run_lock(output), closing(RunStore(output)) as state:
        plan = state.plan
        settings = scoring_settings()
        metadata = plan["metadata"]
        if any(metadata.get(key) != value for key, value in settings.items()) or (
            plan["scorer_fingerprints"] != scorer_fingerprints(metadata["condition_scorers"])
        ):
            raise ValueError("Resume requires the same model, SDK, and scorer settings as the saved run.")
        input_text = (output / "variants.tsv").read_text(encoding="utf-8-sig")
        if hashlib.sha256(input_text.encode("utf-8")).hexdigest() != plan["input_sha256"]:
            raise ValueError("The saved variants.tsv differs from the original run.")
        if cache_dir is None and plan["cache_dir"] is not None:
            cache_dir = Path(plan["cache_dir"])
        return execute(output, state, plan, cache_dir, api_key_file, restore_cache=True)


def main() -> int:
    supported = variant_scorers.get_recommended_scorers(ORGANISM.to_proto())
    choices = [
        name for name, scorer in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items()
        if scorer in supported
    ]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="UTF-8 variant TSV")
    parser.add_argument("--reference", type=Path, help="Indexed GRCh38 FASTA")
    parser.add_argument("--condition", action="append", help="AG2 reference condition; repeatable")
    parser.add_argument("--scorer", action="append", choices=choices, help="AG2 scorer; repeatable")
    parser.add_argument("--panel", type=Path, help="TOML mapping of reference conditions to scorers")
    parser.add_argument("--dry-run", action="store_true", help="Verify REF and count planned score calls offline")
    parser.add_argument("--output-dir", type=Path, help="New directory for results; required for scoring")
    parser.add_argument("--api-key-file", type=Path, help="Private text file containing the API key")
    parser.add_argument("--resume", type=Path, help="Resume a saved scoring run directory")
    parser.add_argument("--cache-dir", type=Path, help="Optional shared AG2 response cache")
    args = parser.parse_args()
    if args.resume is not None:
        if any((args.input, args.reference, args.panel, args.condition, args.scorer, args.output_dir, args.dry_run)):
            parser.error("--resume uses the input and scoring settings saved in the run directory.")
    else:
        if args.input is None or args.reference is None:
            parser.error("--input and --reference are required for a new run.")
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
        if args.resume is not None:
            return resume(args.resume, args.cache_dir, args.api_key_file)
        panel = None
        if args.panel is None:
            condition_scorers = dict.fromkeys(args.condition, args.scorer)
        else:
            name, panel_version, condition_scorers = read_panel(args.panel, choices)
            panel = {"name": name, "version": panel_version, "conditions": condition_scorers}
        return run_batch(args, condition_scorers, panel)
    except Timeout:
        print("The run or cache is locked by another process.", file=sys.stderr)
        return 1
    except (OSError, ValueError, csv.Error, sqlite3.Error) as exc:
        print(exc, file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
