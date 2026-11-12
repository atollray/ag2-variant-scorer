"""Summarize candidate genes from a saved AG2 score.py run without calling the API."""

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import io
import json
import math
from pathlib import Path
import sys


EXPRESSION_SCORERS = ("SC_GENE_EXPRESSION", "RNA_SEQ")
REQUEST_KEY = ("input_id", "requested_condition", "requested_scorer")
INPUT_KEY = ("input_line", "input_id", "requested_scorer")
GENE_COLUMNS = (
    *INPUT_KEY, "gene_id", "gene_name", "rank", "min_raw_score", "max_raw_score",
    "max_abs_raw_score", "meets_threshold", "min_condition", "min_track", "min_imputation_mode",
    "max_condition", "max_track", "max_imputation_mode",
)
CONDITION_STATUSES = (
    "complete", "partial", "api_error", "not_attempted", "no_rows", "no_eligible", "no_finite_scores",
)
VARIANT_COLUMNS = (
    *INPUT_KEY, "evaluation_status", "threshold", "threshold_result", "candidate_genes",
    "conditions_expected", *(f"conditions_{status}" for status in CONDITION_STATUSES),
    "condition_statuses", "rejection_reason",
)


@dataclass(frozen=True)
class InputRow:
    line: int
    input_id: str
    rejection: str | None = None


@dataclass(frozen=True)
class TrackScore:
    value: float
    condition: str
    track: str
    imputation_mode: str
    gene_name: str


def tsv_rows(path: Path, required: tuple[str, ...]):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        headers = reader.fieldnames or []
        if set(required).difference(headers) or len(headers) != len(set(headers)):
            raise ValueError(f"{path.name} requires unique headers including {', '.join(required)}.")
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f"{path.name} line {reader.line_num} does not match the header.")
            yield row


def read_inputs(run_dir: Path) -> list[InputRow]:
    rejected = {}
    rejected_path = run_dir / "rejected.tsv"
    if rejected_path.exists():
        for row in tsv_rows(rejected_path, ("line", "input_id", "error")):
            line = int(row["line"])
            if line in rejected or not row["error"]:
                raise ValueError(f"Invalid rejection for input line {line}.")
            rejected[line] = (row["input_id"], row["error"])

    text = (run_dir / "variants.tsv").read_text(encoding="utf-8-sig")
    # Match the scoring reader's line numbers and its handling of long SV alleles.
    csv.field_size_limit(max(csv.field_size_limit(), len(text)))
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
    headers = reader.fieldnames or []
    if headers.count("input_id") != 1:
        raise ValueError("variants.tsv requires a unique input_id column.")

    inputs = []
    accepted = set()
    for row in reader:
        line = reader.line_num
        input_id = (row["input_id"] or "").strip()
        if line in rejected:
            rejected_id, reason = rejected.pop(line)
            if rejected_id != input_id:
                raise ValueError(f"rejected.tsv input_id does not match variants.tsv line {line}.")
            inputs.append(InputRow(line, input_id, reason))
        else:
            if not input_id or input_id in accepted:
                raise ValueError(f"Missing or duplicate input_id at unrejected input line {line}.")
            accepted.add(input_id)
            inputs.append(InputRow(line, input_id))
    if rejected or not inputs:
        raise ValueError("The input and rejection records do not describe a valid set of input rows.")
    return inputs


def planned_conditions(manifest: dict) -> dict[str, list[str]]:
    plan = manifest.get("condition_scorers")
    if plan is None:
        conditions = manifest.get("conditions")
        scorers = manifest.get("scorers")
        if not isinstance(conditions, list) or not isinstance(scorers, list):
            raise ValueError("request.json is missing its requested conditions or scorers.")
        plan = dict.fromkeys(conditions, scorers)
    if not isinstance(plan, dict):
        raise ValueError("request.json condition_scorers must map conditions to scorer lists.")

    by_scorer = defaultdict(list)
    for condition, scorers in plan.items():
        if not isinstance(scorers, list):
            raise ValueError(f"request.json must list the scorers for {condition}.")
        for name in EXPRESSION_SCORERS:
            if name in scorers:
                by_scorer[name].append(condition)
    if not by_scorer:
        raise ValueError("The run did not request SC_GENE_EXPRESSION or RNA_SEQ.")
    return dict(by_scorer)


def read_requests(run_dir: Path, inputs: list[InputRow], conditions: dict[str, list[str]]):
    accepted = {row.input_id for row in inputs if row.rejection is None}
    expected = {scorer: set(values) for scorer, values in conditions.items()}
    requests = {}
    for row in tsv_rows(run_dir / "requests.tsv", (*REQUEST_KEY, "status", "score_rows")):
        input_id, condition, scorer = key = tuple(row[name] for name in REQUEST_KEY)
        if scorer not in EXPRESSION_SCORERS:
            continue
        if input_id not in accepted or condition not in expected.get(scorer, ()):
            raise ValueError(f"Unplanned expression request: {key}.")
        if key in requests:
            raise ValueError(f"Duplicate expression request: {key}.")
        status = row["status"]
        count = int(row["score_rows"])
        if status not in ("scored", "no_rows", "api_error") or count < 0 or (status == "scored") != (count > 0):
            raise ValueError(f"Invalid request status or score row count: {key}.")
        requests[key] = (status, count)
    return requests


def read_scores(run_dir: Path, requests: dict):
    row_counts = Counter()
    observed = Counter()
    finite = Counter()
    genes = defaultdict(dict)
    fields = (*REQUEST_KEY, "raw_score", "imputation_mode")
    for row in tsv_rows(run_dir / "scores.tsv", fields):
        input_id, condition, scorer = key = tuple(row[name] for name in REQUEST_KEY)
        if scorer not in EXPRESSION_SCORERS:
            continue
        if key not in requests or requests[key][0] != "scored":
            raise ValueError(f"Score row has no matching scored request: {key}.")
        if "gene_id" not in row or "track_name" not in row:
            raise ValueError("scores.tsv expression rows require gene_id and track_name columns.")
        row_counts[key] += 1
        gene_id = row["gene_id"]
        if not gene_id or row["imputation_mode"] != "observed":
            continue
        observed[key] += 1
        value = float(row["raw_score"] or "nan")
        if not math.isfinite(value):
            continue
        finite[key] += 1
        score = TrackScore(value, condition, row["track_name"], row["imputation_mode"], row.get("gene_name") or "")
        group = genes[input_id, scorer]
        low, high = group.get(gene_id, (score, score))
        group[gene_id] = (score if value < low.value else low, score if value > high.value else high)

    for key, (_, expected) in requests.items():
        if row_counts[key] != expected:
            raise ValueError(f"Score row count differs from requests.tsv for {key}: {row_counts[key]} vs {expected}.")
    return genes, observed, finite


def condition_status(request: tuple[str, int], observed_count: int, finite_count: int) -> str:
    status, _ = request
    if status != "scored":
        return status
    if observed_count == 0:
        return "no_eligible"
    if finite_count == 0:
        return "no_finite_scores"
    return "complete" if finite_count == observed_count else "partial"


def write_summaries(output: Path, inputs, conditions, thresholds, requests, genes, observed, finite) -> None:
    with (output / "gene_summary.tsv").open("w", encoding="utf-8", newline="") as gene_handle, (
        output / "variant_summary.tsv"
    ).open("w", encoding="utf-8", newline="") as variant_handle:
        gene_writer = csv.writer(gene_handle, delimiter="\t", lineterminator="\n")
        variant_writer = csv.writer(variant_handle, delimiter="\t", lineterminator="\n")
        gene_writer.writerow(GENE_COLUMNS)
        variant_writer.writerow(VARIANT_COLUMNS)
        for row in inputs:
            for scorer, expected_conditions in conditions.items():
                key = (row.line, row.input_id, scorer)
                threshold = thresholds.get(scorer)
                if row.rejection is not None:
                    variant_writer.writerow((
                        *key, "input_rejected", threshold, "not_evaluated", None,
                        0, *(0 for _ in CONDITION_STATUSES), "{}", row.rejection,
                    ))
                    continue

                statuses = {}
                for condition in expected_conditions:
                    request_key = (row.input_id, condition, scorer)
                    statuses[condition] = condition_status(
                        requests.get(request_key, ("not_attempted", 0)),
                        observed[request_key], finite[request_key],
                    )
                counts = Counter(statuses.values())
                evaluation = "unavailable"
                if counts["complete"] == len(statuses):
                    evaluation = "complete"
                elif counts["complete"] or counts["partial"]:
                    evaluation = "partial"

                effects = []
                for gene_id, (low, high) in genes.get((row.input_id, scorer), {}).items():
                    strongest = low if abs(low.value) >= abs(high.value) else high
                    effects.append((abs(strongest.value), gene_id, strongest.gene_name, low, high))
                effects.sort(key=lambda item: (-item[0], item[1]))
                hits = 0
                for rank, (magnitude, gene_id, gene_name, low, high) in enumerate(effects, start=1):
                    meets_threshold = threshold is not None and magnitude >= threshold
                    gene_writer.writerow((
                        *key, gene_id, gene_name, rank, low.value, high.value, magnitude,
                        str(meets_threshold).lower() if threshold is not None else None,
                        low.condition, low.track, low.imputation_mode,
                        high.condition, high.track, high.imputation_mode,
                    ))
                    if meets_threshold:
                        hits += 1

                result = "not_set"
                if threshold is not None:
                    result = "not_evaluated"
                    if hits:
                        result = "hits"
                    elif evaluation == "complete":
                        result = "no_hits"
                variant_writer.writerow((
                    *key, evaluation, threshold, result, hits if result in ("hits", "no_hits") else None,
                    len(statuses), *(counts[status] for status in CONDITION_STATUSES),
                    json.dumps(statuses, ensure_ascii=False), "",
                ))


def summarize(run_dir: Path, output: Path, thresholds: dict[str, float]) -> None:
    if output.exists():
        raise ValueError(f"Choose a new output directory: {output}")
    manifest = json.loads((run_dir / "request.json").read_text(encoding="utf-8-sig"))
    if not isinstance(manifest, dict):
        raise ValueError("request.json must contain a run record.")
    conditions = planned_conditions(manifest)
    unknown = thresholds.keys() - conditions.keys()
    if unknown:
        raise ValueError(f"Thresholds provided for unrequested scorers: {', '.join(sorted(unknown))}.")
    inputs = read_inputs(run_dir)
    requests = read_requests(run_dir, inputs, conditions)
    genes, observed, finite = read_scores(run_dir, requests)

    output.mkdir(parents=True)
    write_summaries(output, inputs, conditions, thresholds, requests, genes, observed, finite)
    metadata = {
        "source_run": str(run_dir.resolve()),
        "source_request": manifest,
        "expression_conditions_by_scorer": conditions,
        "imputation_mode": "observed",
        "thresholds": thresholds,
    }
    (output / "summary.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True, help="Saved score.py run directory")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory for summaries")
    parser.add_argument(
        "--threshold", action="append", default=[], metavar="SCORER=VALUE",
        help="Nonnegative raw-score magnitude threshold; repeat for different expression scorers",
    )
    args = parser.parse_args()
    thresholds = {}
    try:
        for arg in args.threshold:
            scorer, _, value = arg.partition("=")
            if scorer not in EXPRESSION_SCORERS or scorer in thresholds:
                raise ValueError(f"Unknown or duplicate threshold scorer: {scorer}.")
            threshold = float(value)
            if not math.isfinite(threshold) or threshold < 0:
                raise ValueError(f"Threshold must be finite and nonnegative: {arg}.")
            thresholds[scorer] = threshold
        summarize(args.run_dir, args.output_dir, thresholds)
    except (OSError, ValueError, csv.Error) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"Summaries saved to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
