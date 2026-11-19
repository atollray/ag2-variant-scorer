"""Offline review of gene-linked expression and splicing evidence."""

from collections import Counter, defaultdict
from contextlib import ExitStack, closing, nullcontext
import csv
from dataclasses import asdict, dataclass, fields
from itertools import groupby
import json
import math
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
import tomllib

from .splicing import SPLICE_OUTPUTS, splice_summary
from .storage import run_lock
from .summarize import EXPRESSION_SCORERS, REQUEST_KEY, condition_plan, condition_status, read_inputs, tsv_rows


SC, RNA = EXPRESSION_SCORERS
EXPRESSION = {SC, RNA, SC + "_ACTIVE", RNA + "_ACTIVE"}
KEY = ("input_id", "requested_condition", "scored_interval", "gene_id")
TRACK = ("track_name", "gene_strand", "track_strand", "imputation_mode")
DETAIL = ("gene_name", "requested_scorer", *TRACK, "splice_site_type", "junction_Start", "junction_End",
          "raw_score", "source_record")
MODES = ("observed", "imputed_supported", "imputed_exploratory")
CPM = ("ref_cpm", "alt_cpm", "delta_cpm", "max_cpm", "log1p_max_cpm")
RANK_FIELDS = {"sc": "sc_lfc", "rna": "rna_lfc", "splicing": "splicing_observed",
               "exploratory_splicing": "splicing_imputed_exploratory"}


@dataclass(frozen=True)
class ReviewPolicy:
    min_sc_abs_lfc: float | None = None
    min_sc_active: float | None = None
    rescue_sc_abs_lfc: float | None = None
    min_rna_abs_lfc: float | None = None
    min_splicing: float | None = None
    min_exploratory_splicing: float | None = None
    limit: int | None = None
    conditions: tuple[str, ...] = ()

    @classmethod
    def load(cls, path, plan):
        values = tomllib.loads(path.read_text(encoding="utf-8-sig")) if path else {}
        unknown = values.keys() - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f"Unknown review policy settings: {sorted(unknown)}")
        for key, value in values.items():
            if key == "conditions":
                if not isinstance(value, list) or not all(isinstance(c, str) and c in plan for c in value):
                    raise ValueError("conditions must list conditions in this run")
            elif key == "limit":
                if type(value) is not int or value < 1:
                    raise ValueError("limit must be a positive integer")
            elif type(value) not in (int, float) or not math.isfinite(value) or (key != "min_sc_active" and value < 0):
                raise ValueError(f"Invalid review policy value for {key}")
        policy = cls(**values)
        if policy.min_sc_active is not None and policy.min_sc_abs_lfc is None:
            raise ValueError("min_sc_active requires min_sc_abs_lfc")
        if policy.rescue_sc_abs_lfc is not None and policy.min_sc_active is None:
            raise ValueError("rescue_sc_abs_lfc requires min_sc_active")
        return policy

    def expression_reason(self, scorer, lfc, active):
        if lfc is None:
            return ""
        magnitude = abs(lfc)
        if scorer == RNA:
            return "rna" if self.min_rna_abs_lfc is not None and magnitude >= self.min_rna_abs_lfc else ""
        if self.min_sc_active is not None:
            if active is None:
                return ""
            if active <= self.min_sc_active:
                if self.rescue_sc_abs_lfc is not None and magnitude >= self.rescue_sc_abs_lfc:
                    return "sc_rescue"
                return ""
        return "sc" if self.min_sc_abs_lfc is not None and magnitude >= self.min_sc_abs_lfc else ""


def writer(stack, directory, name, columns):
    handle = stack.enter_context((directory / name).open("w", encoding="utf-8", newline=""))
    result = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    result.writeheader()
    return result


def stage_scores(db, source, plan, inputs, policy, coverage, regions):
    accepted = {item.input_id for item in inputs if item.rejection is None}
    requests = {}
    for r in tsv_rows(source / "requests.tsv", (*REQUEST_KEY, "status", "score_rows", "scored_interval")):
        key = tuple(r[k] for k in REQUEST_KEY)
        vid, cond, scorer = key
        count = int(r["score_rows"])
        if (vid not in accepted or scorer not in plan.get(cond, ()) or key in requests or count < 0
                or r["status"] not in ("scored", "no_rows", "api_error") or (r["status"] == "scored") != (count > 0)):
            raise ValueError(f"Invalid request record for {key}")
        requests[key] = (r["status"], count, r["scored_interval"])
    counts = Counter()
    gene_counts, observed_counts, finite_counts = Counter(), Counter(), Counter()
    columns = (*KEY, *DETAIL)
    db.execute("CREATE TABLE evidence (" + ",".join(columns) + ")")
    insert = "INSERT INTO evidence VALUES (" + ",".join("?" for _ in columns) + ")"
    for record, r in enumerate(tsv_rows(source / "scores.tsv", (*REQUEST_KEY, "raw_score", "imputation_mode")), 1):
        key = tuple(r[k] for k in REQUEST_KEY)
        if key not in requests or requests[key][0] != "scored":
            raise ValueError(f"Unexpected score record {record}")
        counts[key] += 1
        window = requests[key][2]
        if r.get("scored_interval", window) != window:
            raise ValueError(f"Score/request window mismatch at record {record}")
        if policy.conditions and r["requested_condition"] not in policy.conditions:
            continue
        r.update(scored_interval=window, source_record=record)
        scorer = r["requested_scorer"]
        if (scorer not in EXPRESSION and scorer not in SPLICE_OUTPUTS) or not r.get("gene_id"):
            regions.writerow(r)
            continue
        value = float(r["raw_score"] or "nan")
        r["raw_score"] = value if math.isfinite(value) else None
        gene_counts[key] += 1
        if r["imputation_mode"] == "observed":
            observed_counts[key] += 1
            finite_counts[key] += int(r["raw_score"] is not None)
        db.execute(insert, [r.get(k, "") for k in columns])
    if any(counts[k] != entry[1] for k, entry in requests.items()):
        raise ValueError("Score counts differ from the request ledger")
    for item in inputs:
        for condition in policy.conditions or plan:
            for scorer in plan[condition]:
                key = (item.input_id, condition, scorer)
                status, n, window = requests.get(key, ("not_attempted", 0, ""))
                if item.rejection is not None:
                    status, n, window = "input_rejected", 0, ""
                gene_n, obs_n, finite_n = (0, 0, 0) if item.rejection else (
                    gene_counts[key], observed_counts[key], finite_counts[key])
                observed_status = condition_status((status, n), gene_n, obs_n, finite_n)
                coverage.writerow(dict(input_line=item.line, input_id=item.input_id, requested_condition=condition,
                                       requested_scorer=scorer, status=status, score_rows=n, scored_interval=window,
                                       observed_gene_status=observed_status, finite_observed_gene_rows=finite_n))
    return requests


def cpm_values(lfc, active):
    high = math.exp(active) - .001
    low = math.exp(active - abs(lfc)) - .001
    # Float32 log scores can round ln(0.001) slightly down at zero expression.
    ref, alt = (low, high) if lfc >= 0 else (high, low)
    ref, alt = max(0., ref), max(0., alt)
    return dict(zip(CPM, (ref, alt, alt - ref, max(ref, alt), math.log1p(max(ref, alt)))))


def expression_rows(rows, gene, requests, planned, sdk, policy, tracks):
    pairs = defaultdict(dict)
    for r in rows:
        scorer = r["requested_scorer"]
        if scorer not in EXPRESSION:
            continue
        key = (scorer.removesuffix("_ACTIVE"), *(r[k] for k in TRACK))
        if scorer in pairs[key]:
            raise ValueError(f"Ambiguous expression track for {gene['input_id']}, {gene['gene_id']}: {key}")
        pairs[key][scorer] = r
    hits = {}
    for key, pair in sorted(pairs.items()):
        scorer, track, strand, track_strand, mode = key
        active_name = scorer + "_ACTIVE"
        effect_row, active_row = pair.get(scorer), pair.get(active_name)
        lfc = effect_row["raw_score"] if effect_row else None
        active = active_row["raw_score"] if active_row else None
        status = "paired"
        if effect_row is None:
            status = "no_effect_row"
        elif active_row is None:
            request = requests.get((gene["input_id"], gene["requested_condition"], active_name))
            if request:
                status = "no_active_row" if request[0] == "scored" else request[0]
            else:
                status = "not_attempted" if active_name in planned else "active_not_requested"
        elif not all((strand, track_strand, track, mode)):
            status = "missing_track_identity"
            active = None
        elif lfc is None or active is None:
            status = "nonfinite_pair"
        detail = dict(gene, scorer=scorer, track_name=track, gene_strand=strand, track_strand=track_strand,
                      imputation_mode=mode, lfc=lfc, active=active, pair_status=status,
                      effect_record=effect_row["source_record"] if effect_row else None,
                      active_record=active_row["source_record"] if active_row else None)
        if scorer == SC and status == "paired" and sdk == "2.0.0":
            detail.update(cpm_values(lfc, active))
        reason = policy.expression_reason(scorer, lfc, active) if mode == "observed" else ""
        detail["review_reason"] = reason
        tracks.writerow(detail)
        if mode != "observed" or lfc is None:
            continue
        prefix = "sc" if scorer == SC else "rna"
        if reason:
            hits[reason] = max(hits.get(reason, 0.), abs(lfc))
        previous = gene.get(prefix + "_lfc")
        if previous is None or abs(lfc) > abs(previous):
            gene.update({prefix + "_lfc": lfc, prefix + "_track": track, prefix + "_record": detail["effect_record"]})
            if scorer == SC:
                gene.update(sc_active=active, sc_pair_status=status, **{k: detail.get(k) for k in CPM})
    return hits


def splice_rows(rows, gene, components):
    maxima = {}
    for r in rows:
        name, mode = r["requested_scorer"], r["imputation_mode"]
        if name not in SPLICE_OUTPUTS or r["raw_score"] is None:
            continue
        prev = maxima.get((mode, name))
        if prev is None or abs(r["raw_score"]) > abs(prev["raw_score"]):
            maxima[mode, name] = r
    for mode in MODES:
        values = {}
        for name in SPLICE_OUTPUTS:
            r = maxima.get((mode, name))
            if r is not None:
                components.writerow(dict(r))
                values[name] = abs(r["raw_score"])
        summary = splice_summary(values)
        gene["splicing_" + mode] = summary["merged_splicing"]
        if mode == "observed":
            gene.update(values)
            gene["missing_observed_components"] = summary["missing_observed_components"]


def assign_ranks(genes):
    for channel, column in RANK_FIELDS.items():
        values = sorted({abs(g[column]) for g in genes if g.get(column) is not None}, reverse=True)
        ranks = {v: i for i, v in enumerate(values, 1)}
        for g in genes:
            if g.get(column) is not None:
                g[channel + "_rank"] = ranks[abs(g[column])]


def write_review(pairs, policy, table):
    # Rank each evidence type independently. Different raw-score units never compete directly.
    best_rank = {}
    for channel in (*RANK_FIELDS, "sc_rescue"):
        eligible = [k for k, reasons in pairs.items() if channel in reasons]
        ordered = sorted(eligible, key=lambda k: (-pairs[k][channel], k))
        for rank, key in enumerate(ordered, 1):
            best_rank[key] = min(best_rank.get(key, rank), rank)
    selected = Counter()
    for rank, key in enumerate(sorted(pairs, key=lambda k: (best_rank[k], k)), 1):
        vid, gid, kind = key
        keep = policy.limit is None or rank <= policy.limit
        if keep:
            selected[kind] += 1
        table.writerow(dict(input_id=vid, gene_id=gid, allele_type=kind, review_order=rank,
                            review_reasons=",".join(sorted(pairs[key])), within_limit=int(keep)))
    return selected


def build(source, dest, scratch, policy_path):
    manifest = json.loads((source / "request.json").read_text(encoding="utf-8-sig"))
    plan = condition_plan(manifest)
    policy = ReviewPolicy.load(policy_path, plan)
    inputs = read_inputs(source)
    kinds = {item.input_id: item.allele_type for item in inputs if item.rejection is None}
    pairs = defaultdict(dict)
    with ExitStack() as stack:
        db = stack.enter_context(closing(sqlite3.connect(scratch)))
        db.row_factory = sqlite3.Row
        cover = writer(stack, dest, "coverage.tsv", ("input_line", *REQUEST_KEY, "scored_interval", "status", "score_rows",
                                                    "observed_gene_status", "finite_observed_gene_rows"))
        regions = writer(stack, dest, "other_scores.tsv", (*KEY, *DETAIL))
        requests = stage_scores(db, source, plan, inputs, policy, cover, regions)
        tracks = writer(stack, dest, "expression_tracks.tsv", (*KEY, "scorer", *TRACK, "lfc", "active", *CPM,
                       "pair_status", "effect_record", "active_record", "review_reason"))
        components = writer(stack, dest, "splice_components.tsv", (*KEY, *DETAIL))
        genes_out = writer(stack, dest, "gene_conditions.tsv", (*KEY, "gene_name", "allele_type", "sc_lfc", "sc_track",
                           "sc_record", "sc_active", "sc_pair_status", *CPM, "rna_lfc", "rna_track", "rna_record",
                           *SPLICE_OUTPUTS, *("splicing_" + m for m in MODES), "missing_observed_components",
                           *(c + "_rank" for c in RANK_FIELDS), "review_reasons"))
        review_out = writer(stack, dest, "review_pairs.tsv", ("input_id", "gene_id", "allele_type", "review_order",
                            "review_reasons", "within_limit"))
        rows = stack.enter_context(closing(db.execute(
            "SELECT * FROM evidence ORDER BY input_id, requested_condition, scored_interval, gene_id, source_record"
        )))
        for scope, subset in groupby(rows, key=lambda r: tuple(r[k] for k in KEY[:3])):
            genes = []
            for gid, detail in groupby(subset, key=lambda r: r["gene_id"]):
                group = list(detail)
                gene = dict(zip(KEY, (*scope, gid)), gene_name=group[0]["gene_name"], allele_type=kinds[scope[0]])
                hits = expression_rows(group, gene, requests, plan[scope[1]], manifest.get("sdk_version"), policy, tracks)
                splice_rows(group, gene, components)
                merged = gene["splicing_observed"]
                if policy.min_splicing is not None and merged is not None and merged >= policy.min_splicing:
                    hits["splicing"] = merged
                exploratory = gene["splicing_imputed_exploratory"]
                if (policy.min_exploratory_splicing is not None and exploratory is not None
                        and exploratory >= policy.min_exploratory_splicing):
                    hits["exploratory_splicing"] = exploratory
                gene["review_reasons"] = ",".join(sorted(hits))
                for reason, value in hits.items():
                    reasons = pairs[scope[0], gid, gene["allele_type"]]
                    reasons[reason] = max(reasons.get(reason, 0.), value)
                genes.append(gene)
            assign_ranks(genes)
            genes_out.writerows(genes)
        selected = write_review(pairs, policy, review_out)
    info = dict(source_run=str(source.resolve()), source_request=manifest, review_policy=asdict(policy),
                eligible_unique_pairs=len(pairs), selected_unique_pairs=sum(selected.values()),
                selected_pairs_by_allele_type=selected)
    (dest / "summary.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")


def write_candidates(source: Path, output: Path, policy_path: Path | None):
    if output.exists():
        raise ValueError(f"Choose a new output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as tmp, (
        run_lock(source) if (source / "run.sqlite").is_file() else nullcontext()
    ):
        stage = Path(tmp) / "summary"
        stage.mkdir()
        build(source, stage, Path(tmp) / "rows.sqlite", policy_path)
        stage.rename(output)
