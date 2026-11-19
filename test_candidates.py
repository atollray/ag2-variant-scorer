"""Offline CLI regression cases with synthetic scores."""

import csv
import io
import json
import math
import tempfile
import unittest
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from ag2_variant_scorer.cli import main
from ag2_variant_scorer.splicing import SPLICE_OUTPUTS

SC = "SC_GENE_EXPRESSION"
ACTIVE = SC + "_ACTIVE"
SITES, USAGE, JUNCTIONS = SPLICE_OUTPUTS


class CandidateTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.run = self.root / "run"
        self.run.mkdir()
        self.out = self.root / "out"

    def write(self, name, header, rows):
        with (self.run / name).open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)

    def fixture(self, plan, rows, sdk="2.0.0"):
        # rows: condition, scorer, gene, track, raw, mode and optional genomic metadata.
        prepared = []
        for cond, scorer, gene, track, raw, mode, *extra in rows:
            strand, ts, start, end = extra or ("+", ".", "", "")
            prepared.append(("v", cond, scorer, "window", gene, gene, track, strand, ts, raw, mode, start, end))
        counts = Counter(r[:3] for r in prepared)
        self.write("variants.tsv", ["input_id", "ref", "alt"], [("v", "A", "T")])
        self.write("requests.tsv", ["input_id", "requested_condition", "requested_scorer", "scored_interval", "status", "score_rows"],
                   [("v", c, s, "window", "scored" if counts["v", c, s] else "no_rows", counts["v", c, s])
                    for c, scorers in plan.items() for s in scorers])
        self.write("scores.tsv", ["input_id", "requested_condition", "requested_scorer", "scored_interval",
                   "gene_id", "gene_name", "track_name", "gene_strand", "track_strand", "raw_score",
                   "imputation_mode", "junction_Start", "junction_End"], prepared)
        (self.run / "request.json").write_text(json.dumps(dict(condition_scorers=plan, sdk_version=sdk)))

    def invoke(self, settings=None):
        args = ["summarize", "--gene-conditions", "--run-dir", str(self.run), "--output-dir", str(self.out)]
        if settings is not None:
            policy = self.root / "policy.toml"
            policy.write_text(settings)
            args += ["--review-policy", str(policy)]
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            return main(args)

    def read(self, name):
        with (self.out / (name + ".tsv")).open(encoding="utf-8") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_splice_only_per_gene_condition_mode_and_junction(self):
        rows = [
            ("a", SITES, "G1", "donor", .9, "observed"),
            ("a", USAGE, "G1", "usage", .1, "observed"),
            ("a", JUNCTIONS, "G1", "junction", 1, "observed", "+", ".", "30", "50"),
            ("a", JUNCTIONS, "G2", "junction", 5, "observed", "-", ".", "10", "20"),
            ("b", USAGE, "G2", "usage", .8, "observed"),
            ("a", USAGE, "G2", "usage", .7, "imputed_exploratory"),
            ("a", SITES, "G2", "site", .4, "imputed_exploratory"),
            ("a", JUNCTIONS, "G2", "j", .5, "imputed_exploratory"),
            ("a", "CAGE", "G2", "center", 20, "observed"),
        ]
        self.fixture({"a": [*SPLICE_OUTPUTS, "CAGE"], "b": list(SPLICE_OUTPUTS)}, rows)
        self.assertEqual(self.invoke("min_splicing = 1.0\nmin_exploratory_splicing = 1.0"), 0)
        g1, g2, other = self.read("gene_conditions")
        self.assertAlmostEqual(float(g1["splicing_observed"]), 1.2)
        self.assertEqual(g1["review_reasons"], "splicing")
        self.assertEqual(g2["splicing_observed"], "")
        self.assertEqual(g2["missing_observed_components"], SITES + "," + USAGE)
        self.assertEqual(g2["review_reasons"], "exploratory_splicing")
        self.assertEqual(other["splicing_observed"], "")
        j = next(r for r in self.read("splice_components") if r["source_record"] == "4")
        self.assertEqual((j["gene_id"], j["gene_strand"], j["junction_Start"], j["junction_End"]), ("G2", "-", "10", "20"))
        self.assertEqual(self.read("other_scores")[0]["source_record"], "9")

    def test_track_pairing_selection_and_distinct_pair_limit(self):
        rows = [
            ("a", SC, "G", "low", -.4, "observed"),
            ("a", ACTIVE, "G", "low", -2, "observed"),
            ("a", SC, "G", "high", -.2, "observed"),
            ("a", ACTIVE, "G", "high", 2, "observed"),
            ("a", SC, "H", "x", 2, "observed"),
            ("a", ACTIVE, "H", "x", 1, "observed", "-", ".", "", ""),
            ("b", SC, "G", "x", 1, "observed"),
            ("b", ACTIVE, "G", "x", -.1, "observed"),
            ("b", SC, "Z", "x", .15, "observed"),
            ("b", ACTIVE, "Z", "x", 1.5, "observed"),
        ]
        self.fixture({"a": [SC, ACTIVE], "b": [SC, ACTIVE]}, rows)
        self.assertEqual(self.invoke("min_sc_abs_lfc=0.1\nmin_sc_active=1\nrescue_sc_abs_lfc=0.5\nlimit=1"), 0)
        g, h, _, _ = self.read("gene_conditions")
        self.assertEqual(g["review_reasons"], "sc")
        self.assertEqual(float(g["sc_lfc"]), -.4)
        self.assertAlmostEqual(float(g["ref_cpm"]), math.exp(-2)-.001)
        self.assertAlmostEqual(float(g["alt_cpm"]), math.exp(-2.4)-.001)
        self.assertEqual(h["sc_active"], "")
        self.assertEqual(h["review_reasons"], "")
        hits = self.read("review_pairs")
        self.assertEqual([(r["gene_id"], r["within_limit"]) for r in hits], [("G", "1"), ("Z", "0")])
        self.assertEqual(hits[0]["review_reasons"], "sc,sc_rescue")

    def test_no_active_legacy_and_nonfinite(self):
        self.fixture({"a": [SC], "b": [SC, ACTIVE]}, [
            ("a", SC, "G", "x", -.7, "observed"), ("b", SC, "G", "x", -.7, "observed"),
            ("b", ACTIVE, "G", "x", "NaN", "observed")])
        self.assertEqual(self.invoke(), 0)
        a, b = self.read("gene_conditions")
        self.assertEqual((a["sc_lfc"], a["sc_pair_status"], b["sc_pair_status"]), ("-0.7", "active_not_requested", "nonfinite_pair"))
        self.assertEqual((a["ref_cpm"], b["ref_cpm"]), ("", ""))

    def test_unknown_sdk_retains_raw_and_condition_filter(self):
        self.fixture({"a": [SC, ACTIVE], "b": [SC]}, [("a", SC, "G", "x", 1, "observed"),
                     ("a", ACTIVE, "G", "x", 3, "observed"), ("b", SC, "H", "x", 4, "observed")], sdk="unknown")
        self.assertEqual(self.invoke('conditions=["a"]'), 0)
        gene, = self.read("gene_conditions")
        self.assertEqual((gene["sc_active"], gene["ref_cpm"]), ("3.0", ""))

    def test_no_rows_and_rejected_policy(self):
        self.fixture({"a": [SC]}, [])
        self.assertEqual(self.invoke("typo=1"), 1)
        self.assertFalse(self.out.exists())
        self.assertEqual(self.invoke(), 0)
        self.assertEqual(self.read("gene_conditions"), [])
        self.assertEqual(self.read("coverage")[0]["status"], "no_rows")

    def test_window_mismatch_fails_without_partial_output(self):
        self.fixture({"a": [SC]}, [("a", SC, "G", "t", .4, "observed")])
        path = self.run / "scores.tsv"
        path.write_text(path.read_text().replace("\twindow\t", "\tother\t"))
        self.assertEqual(self.invoke(), 1)
        self.assertFalse(self.out.exists())

    def test_rejected_malformed_row_and_multiple_reasons_count_once(self):
        plan = {"a": [SC, *SPLICE_OUTPUTS]}
        self.fixture(plan, [("a", scorer, "G", "track", 1, "observed") for scorer in plan["a"]])
        with (self.run / "variants.tsv").open("a") as f:
            f.write("v\n")
        self.write("rejected.tsv", ["line", "input_id", "error"], [(3, "v", "Missing allele")])
        self.assertEqual(self.invoke("min_sc_abs_lfc=0.1\nmin_splicing=0.1"), 0)
        pair, = self.read("review_pairs")
        self.assertEqual(pair["review_reasons"], "sc,splicing")
        bad = [r for r in self.read("coverage") if r["input_line"] == "3"]
        self.assertEqual(len(bad), 4)
        self.assertTrue(all(r["observed_gene_status"] == "input_rejected" for r in bad))
