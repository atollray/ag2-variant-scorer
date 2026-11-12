"""Offline regression tests for gene summaries from saved scoring runs."""

from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

from ag2_variant_scorer import summarize


CONE = "UBERON:0000053_CL:0000573"
RPE = "UBERON:0000053_CL:0002586"
PERIPHERAL = "UBERON:0013682_CL:0002586"
RNA = "RNA_SEQ"
SC = "SC_GENE_EXPRESSION"
INTERVAL = "chr22:35677410-36725986:."


class SummarizeCliTest(unittest.TestCase):
    def setUp(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.run_dir = directory / "run"
        self.run_dir.mkdir()
        self.output = directory / "summary"

    def write_tsv(self, name, columns, rows):
        with (self.run_dir / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            writer.writerows(rows)

    def write_run(self, input_ids, condition_scorers, requests, scores, *, legacy=False, rejected=()):
        self.write_tsv(
            "variants.tsv",
            ("input_id", "chromosome", "position", "ref", "alt"),
            ((input_id, "chr22", 36201698 + offset, "A", "C")
             for offset, input_id in enumerate(input_ids)),
        )
        score_counts = Counter(row[:3] for row in scores)
        self.write_tsv(
            "requests.tsv",
            ("input_id", "requested_condition", "requested_scorer", "scored_interval", "status",
             "score_rows", "error"),
            ((*key, INTERVAL, status, score_counts[key], "UNAVAILABLE" if status == "api_error" else "")
             for key, status in requests),
        )
        # The scorer writes a shorter header when it has no returned score rows.
        columns = ("input_id", "requested_condition", "requested_scorer")
        if scores:
            columns += ("gene_id", "gene_name", "track_name")
        self.write_tsv("scores.tsv", columns + ("raw_score", "imputation_mode"), scores)
        if rejected:
            self.write_tsv("rejected.tsv", ("line", "input_id", "error"), rejected)

        manifest = {
            "assembly": "GRCh38",
            "model_version": "ALL_FOLDS",
            "sdk_version": "2.0.0",
            "window_bp": 1048576,
            "conditions": list(condition_scorers),
            "scorers": sorted({scorer for scorers in condition_scorers.values() for scorer in scorers}),
            "verified_variants": len(input_ids) - len(rejected),
            "rejected_rows": len(rejected),
            "planned_score_calls": (len(input_ids) - len(rejected)) * len(condition_scorers),
            "attempted_score_calls": len({key[:2] for key, _ in requests}),
        }
        if not legacy:
            manifest["condition_scorers"] = condition_scorers
        (self.run_dir / "request.json").write_text(json.dumps(manifest), encoding="utf-8")

    def run_cli(self, *thresholds):
        argv = ["summarize.py", "--run-dir", str(self.run_dir), "--output-dir", str(self.output)]
        for threshold in thresholds:
            argv.extend(("--threshold", threshold))
        stdout, stderr = io.StringIO(), io.StringIO()
        with mock.patch.object(sys, "argv", argv), redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = summarize.main()
        self.assertEqual(exit_code, 0, stderr.getvalue())

    def read_tsv(self, name):
        with (self.output / name).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_legacy_manifest_preserves_signed_extremes_and_ranks_by_scorer(self):
        self.write_run(
            ["001"],
            {CONE: [RNA, SC], RPE: [RNA, SC]},
            [(("001", condition, scorer), "scored") for condition in (CONE, RPE) for scorer in (RNA, SC)],
            [
                ("001", CONE, SC, "ENSG1", "FIRST", "cone sc", -1.25, "observed"),
                ("001", RPE, SC, "ENSG1", "FIRST", "rpe sc", 0.8, "observed"),
                ("001", RPE, SC, "ENSG1", "FIRST", "exploratory", 20.0, "imputed_exploratory"),
                ("001", RPE, SC, "ENSG2", "SECOND", "rpe sc", 0.95, "observed"),
                ("001", CONE, RNA, "ENSG1", "FIRST", "cone bulk", 0.05, "observed"),
                ("001", RPE, RNA, "ENSG1", "FIRST", "rpe bulk", 0.5, "observed"),
                ("001", RPE, RNA, "ENSG2", "SECOND", "rpe bulk", 1.1, "observed"),
            ],
            legacy=True,
        )

        self.run_cli(f"{RNA}=1.2")

        genes = {(row["requested_scorer"], row["gene_id"]): row for row in self.read_tsv("gene_summary.tsv")}
        self.assertEqual(len(genes), 4)
        sc_first = genes[SC, "ENSG1"]
        self.assertEqual(sc_first["input_id"], "001")
        self.assertEqual(float(sc_first["min_raw_score"]), -1.25)
        self.assertEqual(float(sc_first["max_raw_score"]), 0.8)
        self.assertEqual(float(sc_first["max_abs_raw_score"]), 1.25)
        self.assertEqual((sc_first["min_condition"], sc_first["min_track"], sc_first["min_imputation_mode"]),
                         (CONE, "cone sc", "observed"))
        self.assertEqual((sc_first["max_condition"], sc_first["max_track"], sc_first["max_imputation_mode"]),
                         (RPE, "rpe sc", "observed"))
        self.assertEqual(int(sc_first["rank"]), 1)
        self.assertEqual(int(genes[RNA, "ENSG1"]["rank"]), 2)
        self.assertEqual(sc_first["meets_threshold"], "")
        self.assertEqual(genes[RNA, "ENSG1"]["meets_threshold"], "false")

        variants = {row["requested_scorer"]: row for row in self.read_tsv("variant_summary.tsv")}
        self.assertEqual((variants[SC]["evaluation_status"], variants[SC]["threshold_result"]),
                         ("complete", "not_set"))
        self.assertEqual((variants[RNA]["evaluation_status"], variants[RNA]["threshold_result"]),
                         ("complete", "no_hits"))
        self.assertEqual(int(variants[RNA]["candidate_genes"]), 0)

    def test_panel_reports_missing_conditions_without_dropping_partial_hits_or_valid_duplicate_id(self):
        requests = [
            (("complete", CONE, RNA), "scored"),
            (("complete", CONE, SC), "scored"),
            (("complete", RPE, RNA), "scored"),
            (("complete", PERIPHERAL, RNA), "scored"),
            (("absent", CONE, RNA), "no_rows"),
            (("absent", CONE, SC), "no_rows"),
            (("absent", RPE, RNA), "scored"),
            (("absent", PERIPHERAL, RNA), "scored"),
            (("partial", CONE, RNA), "scored"),
            (("partial", CONE, SC), "scored"),
            (("partial", RPE, RNA), "api_error"),
        ]
        self.write_run(
            ["complete", "complete", "absent", "partial", "not_started"],
            {CONE: [RNA, SC], RPE: [RNA], PERIPHERAL: [RNA]},
            requests,
            [
                ("complete", CONE, RNA, "ENSG1", "FIRST", "bulk A", 0.0, "observed"),
                ("complete", CONE, SC, "ENSG1", "FIRST", "sc", 0.1, "observed"),
                ("complete", RPE, RNA, "ENSG1", "FIRST", "bulk B", -0.1, "observed"),
                ("complete", PERIPHERAL, RNA, "ENSG1", "FIRST", "bulk C", 0.2, "observed"),
                ("absent", RPE, RNA, "ENSG1", "FIRST", "exploratory", -4.0, "imputed_exploratory"),
                ("absent", PERIPHERAL, RNA, "ENSG1", "FIRST", "nonfinite", "inf", "observed"),
                ("absent", PERIPHERAL, RNA, "ENSG1", "FIRST", "missing", "", "observed"),
                ("partial", CONE, RNA, "ENSG1", "FIRST", "bulk A", -0.5, "observed"),
                ("partial", CONE, SC, "ENSG1", "FIRST", "sc", 0.2, "observed"),
            ],
            rejected=[(3, "complete", "Duplicate input_id.")],
        )

        self.run_cli(f"{RNA}=0.5", f"{SC}=0.5")

        variants = {(int(row["input_line"]), row["requested_scorer"]): row
                    for row in self.read_tsv("variant_summary.tsv")}
        self.assertEqual(len(variants), 10)
        self.assertEqual((variants[2, RNA]["evaluation_status"], variants[2, RNA]["threshold_result"]),
                         ("complete", "no_hits"))
        self.assertEqual((variants[3, RNA]["evaluation_status"], variants[3, RNA]["threshold_result"]),
                         ("input_rejected", "not_evaluated"))
        self.assertEqual(variants[3, RNA]["rejection_reason"], "Duplicate input_id.")

        absent = variants[4, RNA]
        self.assertEqual((absent["evaluation_status"], absent["threshold_result"]),
                         ("unavailable", "not_evaluated"))
        self.assertEqual(json.loads(absent["condition_statuses"]),
                         {CONE: "no_rows", RPE: "no_eligible", PERIPHERAL: "no_finite_scores"})

        partial = variants[5, RNA]
        self.assertEqual((partial["evaluation_status"], partial["threshold_result"]), ("partial", "hits"))
        self.assertEqual(int(partial["candidate_genes"]), 1)
        self.assertEqual(json.loads(partial["condition_statuses"]),
                         {CONE: "complete", RPE: "api_error", PERIPHERAL: "not_attempted"})
        self.assertEqual((variants[5, SC]["evaluation_status"], variants[5, SC]["threshold_result"]),
                         ("complete", "no_hits"))
        self.assertEqual(int(variants[5, SC]["conditions_expected"]), 1)

        never_started = variants[6, RNA]
        self.assertEqual((never_started["evaluation_status"], never_started["threshold_result"]),
                         ("unavailable", "not_evaluated"))
        self.assertEqual(int(never_started["conditions_not_attempted"]), 3)

        genes = self.read_tsv("gene_summary.tsv")
        self.assertEqual({int(row["input_line"]) for row in genes}, {2, 5})
        self.assertEqual([(row["input_id"], row["requested_scorer"]) for row in genes
                          if row["meets_threshold"] == "true"], [("partial", RNA)])

    def test_nonfinite_observed_values_prevent_reporting_no_hits_from_partial_scores(self):
        self.write_run(
            ["partial"], {CONE: [RNA]}, [(("partial", CONE, RNA), "scored")],
            [
                ("partial", CONE, RNA, "ENSG1", "FIRST", "finite", 0.1, "observed"),
                ("partial", CONE, RNA, "ENSG2", "SECOND", "nonfinite", "inf", "observed"),
            ],
        )

        self.run_cli(f"{RNA}=0.5")

        row, = self.read_tsv("variant_summary.tsv")
        self.assertEqual((row["evaluation_status"], row["threshold_result"]), ("partial", "not_evaluated"))
        self.assertEqual(json.loads(row["condition_statuses"]), {CONE: "partial"})
        gene, = self.read_tsv("gene_summary.tsv")
        self.assertEqual(gene["gene_id"], "ENSG1")

    def test_no_returned_scores_accepts_empty_header_without_gene_columns(self):
        self.write_run(["empty"], {CONE: [RNA]}, [(("empty", CONE, RNA), "no_rows")], [])

        self.run_cli(f"{RNA}=0.5")

        self.assertEqual(self.read_tsv("gene_summary.tsv"), [])
        row, = self.read_tsv("variant_summary.tsv")
        self.assertEqual((row["evaluation_status"], row["threshold_result"]), ("unavailable", "not_evaluated"))
        self.assertEqual(int(row["conditions_no_rows"]), 1)

    def test_scored_rows_without_gene_ids_are_distinct_from_unobserved_gene_rows(self):
        self.write_run(
            ["v1"], {CONE: [RNA], RPE: [RNA]},
            [(("v1", CONE, RNA), "scored"), (("v1", RPE, RNA), "scored")],
            [
                ("v1", CONE, RNA, "", "", "observed track", 0.75, "observed"),
                ("v1", RPE, RNA, "ENSG1", "FIRST", "imputed track", 0.75, "imputed_supported"),
            ],
        )

        self.run_cli(f"{RNA}=0.5")

        self.assertEqual(self.read_tsv("gene_summary.tsv"), [])
        row, = self.read_tsv("variant_summary.tsv")
        self.assertEqual((row["evaluation_status"], row["threshold_result"]), ("unavailable", "not_evaluated"))
        self.assertEqual(json.loads(row["condition_statuses"]), {CONE: "no_gene_rows", RPE: "no_eligible"})
        self.assertEqual(int(row["conditions_no_gene_rows"]), 1)
        self.assertEqual(int(row["conditions_no_eligible"]), 1)
