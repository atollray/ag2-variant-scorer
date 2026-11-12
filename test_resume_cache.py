"""Exercise the scoring CLI's restart and cache behavior without API calls."""

from contextlib import redirect_stdout
import csv
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

import score
from test_score import CONDITIONS, Unavailable, score_row


RNA = "RNA_SEQ"
ACTIVE = "RNA_SEQ_ACTIVE"
SC = "SC_GENE_EXPRESSION"


class AbruptStop(BaseException):
    """Stop execution outside the CLI's normal error handling."""


class ResumeCacheCliTest(unittest.TestCase):
    def setUp(self):
        self.directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.reference = self.directory / "reference.fa"
        self.variants = self.directory / "variants.tsv"
        self.cache = self.directory / "cache"
        bases = b"A" * 1_200_000
        self.reference.write_bytes(b">chr22\n" + bases + b"\n")
        self.reference.with_suffix(".fa.fai").write_bytes(
            f"chr22\t{len(bases)}\t7\t{len(bases)}\t{len(bases) + 1}\n".encode("ascii")
        )
        self.write_variants("v1\tchr22\t600001\tA\tC\n")
        self.create = self.enterContext(mock.patch.object(score.dna_client, "create"))
        self.model = self.create.return_value
        self.model.supported_conditions.return_value = pd.DataFrame({
            "condition": CONDITIONS,
            "biosample_name": ["macula cone cell", "macula retinal pigment epithelial cell"],
            "rna_seq": ["observed", "imputed_supported"],
        })
        self.enterContext(mock.patch.object(
            score.variant_scorers, "tidy_anndata", side_effect=lambda table, **_: table.copy()
        ))

    def write_variants(self, rows: str):
        header = "input_id\tchromosome\tposition\tref\talt\n"
        self.variants.write_text(header + rows, encoding="utf-8")

    def invoke(self, *options: str) -> int:
        with mock.patch.object(sys, "argv", ["score.py", *options]), mock.patch.dict(
            os.environ, {"ALPHAGENOME_API_KEY": "offline-fixture-key"}
        ):
            return score.main()

    def start(
        self, output: Path | None, conditions: tuple[str, ...], scorers: tuple[str, ...],
        *, use_cache: bool = False, dry_run: bool = False,
    ) -> int:
        options = ["--input", str(self.variants), "--reference", str(self.reference)]
        for condition in conditions:
            options.extend(("--condition", condition))
        for scorer in scorers:
            options.extend(("--scorer", scorer))
        if output is not None:
            options.extend(("--output-dir", str(output)))
        if use_cache:
            options.extend(("--cache-dir", str(self.cache)))
        if dry_run:
            options.append("--dry-run")
        return self.invoke(*options)

    def read_rows(self, directory: Path, name: str) -> list[dict[str, str]]:
        with (directory / name).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def assert_call_counts(self, directory: Path, *, planned: int, attempted: int):
        manifest = json.loads((directory / "request.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["planned_score_calls"], planned)
        self.assertEqual(manifest["attempted_score_calls"], attempted)

    def test_resume_retries_only_failed_and_unattempted_groups(self):
        output = self.directory / "interrupted-run"
        self.write_variants("v1\tchr22\t600001\tA\tC\nv2\tchr22\t600002\tA\tT\n")
        self.model.score_variant.side_effect = [
            [score_row(-0.5, "observed"), pd.DataFrame()],
            Unavailable("Offline fixture RPC failure."),
        ]
        self.assertEqual(self.start(output, CONDITIONS, (RNA, ACTIVE)), 1)
        self.assert_call_counts(output, planned=4, attempted=2)
        self.assertEqual([row["status"] for row in self.read_rows(output, "requests.tsv")],
                         ["scored", "no_rows", "api_error", "api_error"])

        # The saved run is sufficient even when the original local inputs are gone.
        self.variants.unlink()
        self.reference.unlink()
        self.model.reset_mock()
        self.model.score_variant.side_effect = [
            [score_row(value, "observed"), score_row(2.0, "observed")]
            for value in (0.75, 1.0, 1.25)
        ]
        self.assertEqual(self.invoke("--resume", str(output)), 0)
        self.assertEqual([
            (call.kwargs["variant"].position, call.kwargs["condition"])
            for call in self.model.score_variant.call_args_list
        ], [(600001, CONDITIONS[1]), (600002, CONDITIONS[0]), (600002, CONDITIONS[1])])
        self.assert_call_counts(output, planned=4, attempted=5)
        requests = self.read_rows(output, "requests.tsv")
        self.assertEqual([row["status"] for row in requests], ["scored", "no_rows"] + ["scored"] * 6)
        original = self.read_rows(output, "scores.tsv")[0]
        self.assertEqual((original["input_id"], float(original["raw_score"])), ("v1", -0.5))

    def test_new_runs_reuse_each_scorer_and_fetch_missing_scorers_together(self):
        first = CONDITIONS[0]
        original = self.directory / "original"
        self.model.score_variant.return_value = [score_row(-0.75, "observed")]
        self.assertEqual(self.start(original, (first,), (RNA,), use_cache=True), 0)
        self.assert_call_counts(original, planned=1, attempted=1)

        mixed = self.directory / "mixed"
        self.model.reset_mock()
        self.model.score_variant.return_value = [pd.DataFrame(), score_row(1.25, "observed", SC)]
        self.assertEqual(self.start(mixed, (first,), (RNA, ACTIVE, SC), use_cache=True), 0)
        self.model.score_variant.assert_called_once()
        self.assertEqual(list(self.model.score_variant.call_args.kwargs["variant_scorers"]), [
            score.variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in (ACTIVE, SC)
        ])
        self.assert_call_counts(mixed, planned=1, attempted=1)
        mixed_requests = {row["requested_scorer"]: row for row in self.read_rows(mixed, "requests.tsv")}
        self.assertEqual(mixed_requests[RNA]["source"], "cache")
        self.assertEqual(mixed_requests[ACTIVE]["status"], "no_rows")
        self.assertEqual(mixed_requests[SC]["source"], "api")

        self.write_variants("renamed\tchr22\t600001\tA\tC\n")
        cached = self.directory / "cached"
        panel = self.directory / "reuse.toml"
        panel.write_text(
            f'name = "reuse"\nversion = 1\n[conditions]\n{json.dumps(first)} = {json.dumps([SC, RNA])}\n',
            encoding="utf-8",
        )
        self.create.reset_mock()
        with mock.patch.object(score, "api_key_from_file_or_environment") as key:
            self.assertEqual(self.invoke(
                "--input", str(self.variants), "--reference", str(self.reference),
                "--panel", str(panel), "--cache-dir", str(self.cache), "--output-dir", str(cached),
            ), 0)
        key.assert_not_called()
        self.create.assert_not_called()
        self.assert_call_counts(cached, planned=1, attempted=0)
        self.assertEqual({row["source"] for row in self.read_rows(cached, "requests.tsv")}, {"cache"})
        self.assertEqual({
            (row["input_id"], row["requested_scorer"]): float(row["raw_score"])
            for row in self.read_rows(cached, "scores.tsv")
        }, {("renamed", RNA): -0.75, ("renamed", SC): 1.25})
        condition, = self.read_rows(cached, "conditions.tsv")
        self.assertEqual((condition["condition"], condition["rna_seq"]), (first, "observed"))

        # Empty success is complete within a run, but is queried again in a new run.
        queried = self.directory / "query-empty-again"
        self.model.score_variant.return_value = [score_row(5.0, "observed")]
        self.assertEqual(self.start(queried, (first,), (ACTIVE,), use_cache=True), 0)
        self.model.score_variant.assert_called_once()
        self.assert_call_counts(queried, planned=1, attempted=1)
        self.assertEqual(self.read_rows(queried, "requests.tsv")[0]["source"], "api")

    def test_cache_aware_dry_run_is_offline_and_reports_only_uncached_calls(self):
        first, second = CONDITIONS
        stdout = io.StringIO()
        with redirect_stdout(stdout), mock.patch.object(score, "api_key_from_file_or_environment") as key:
            self.assertEqual(self.start(None, (first,), (RNA,), use_cache=True, dry_run=True), 0)
        key.assert_not_called()
        self.create.assert_not_called()
        self.assertFalse(self.cache.exists())
        self.assertIn("Uncached score calls (upper bound): 1", stdout.getvalue().splitlines())

        seed = self.directory / "seed"
        self.model.score_variant.return_value = [score_row(0.4, "observed")]
        self.assertEqual(self.start(seed, (first,), (RNA,), use_cache=True), 0)
        self.create.reset_mock()
        self.write_variants("v1\tchr22\t600001\tA\tC\ndifferent\tchr22\t600001\tA\tT\n")
        stdout = io.StringIO()
        with redirect_stdout(stdout), mock.patch.object(score, "api_key_from_file_or_environment") as key:
            self.assertEqual(self.start(None, (first, second), (RNA,), use_cache=True, dry_run=True), 0)
        key.assert_not_called()
        self.create.assert_not_called()
        self.assertIn("Planned score calls: 4", stdout.getvalue().splitlines())
        self.assertIn("Uncached score calls (upper bound): 3", stdout.getvalue().splitlines())

    def test_cache_aware_dry_run_bounds_calls_for_duplicate_alleles(self):
        condition = CONDITIONS[0]
        self.write_variants("v1\tchr22\t600001\tA\tC\nv2\tchr22\t600001\tA\tC\n")
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            self.assertEqual(self.start(None, (condition,), (RNA,), use_cache=True, dry_run=True), 0)
        self.create.assert_not_called()
        self.assertIn("Uncached score calls (upper bound): 2", stdout.getvalue().splitlines())

        output = self.directory / "duplicate-alleles"
        self.model.score_variant.return_value = [score_row(0.125, "observed")]
        self.assertEqual(self.start(output, (condition,), (RNA,), use_cache=True), 0)
        self.model.score_variant.assert_called_once()
        self.assert_call_counts(output, planned=2, attempted=1)
        self.assertEqual({
            row["input_id"]: row["source"] for row in self.read_rows(output, "requests.tsv")
        }, {"v1": "api", "v2": "cache"})

    def test_resume_repopulates_cache_after_a_stopped_shared_write(self):
        condition = CONDITIONS[0]
        output = self.directory / "interrupted-cache-write"
        self.model.score_variant.return_value = [score_row(0.25, "observed")]
        with mock.patch.object(score.SharedCache, "save", side_effect=AbruptStop()) as save:
            with self.assertRaises(AbruptStop):
                self.start(output, (condition,), (RNA,), use_cache=True)
        self.model.score_variant.assert_called_once()
        save.assert_called_once()

        peer = self.directory / "peer-after-resume"
        self.create.reset_mock()
        with mock.patch.object(score, "api_key_from_file_or_environment") as key:
            self.assertEqual(self.invoke("--resume", str(output)), 0)
            self.assertEqual(self.start(peer, (condition,), (RNA,), use_cache=True), 0)
        key.assert_not_called()
        self.create.assert_not_called()
        self.assert_call_counts(output, planned=1, attempted=1)
        self.assert_call_counts(peer, planned=1, attempted=0)
        self.assertEqual(self.read_rows(output, "requests.tsv")[0]["source"], "api")
        self.assertEqual(self.read_rows(peer, "requests.tsv")[0]["source"], "cache")
        self.assertEqual(self.read_rows(peer, "scores.tsv"), self.read_rows(output, "scores.tsv"))

    def test_resume_uses_its_saved_cache_dir_for_a_failed_condition(self):
        first, second = CONDITIONS
        output = self.directory / "interrupted-cache-run"
        self.model.score_variant.side_effect = [
            [score_row(0.125, "observed")], Unavailable("Offline fixture RPC failure.")
        ]
        self.assertEqual(self.start(output, CONDITIONS, (RNA,), use_cache=True), 1)

        self.model.score_variant.side_effect = [[score_row(0.25, "observed")]]
        peer = self.directory / "peer"
        self.assertEqual(self.start(peer, (second,), (RNA,), use_cache=True), 0)

        self.create.reset_mock()
        with mock.patch.object(score, "api_key_from_file_or_environment") as key:
            self.assertEqual(self.invoke("--resume", str(output)), 0)
        key.assert_not_called()
        self.create.assert_not_called()
        self.assert_call_counts(output, planned=2, attempted=2)
        self.assertEqual({
            row["requested_condition"]: (row["status"], row["source"])
            for row in self.read_rows(output, "requests.tsv")
        }, {first: ("scored", "api"), second: ("scored", "cache")})

    def test_resume_recovers_a_saved_response_after_an_abrupt_stop(self):
        first, second = CONDITIONS
        output = self.directory / "abrupt-run"
        self.model.score_variant.side_effect = [[score_row(-0.125, "observed")], AbruptStop()]
        with self.assertRaises(AbruptStop):
            self.start(output, CONDITIONS, (RNA,))

        self.model.reset_mock()
        self.model.score_variant.side_effect = [[score_row(0.25, "observed")]]
        self.assertEqual(self.invoke("--resume", str(output)), 0)
        self.model.score_variant.assert_called_once()
        self.assertEqual(self.model.score_variant.call_args.kwargs["condition"], second)
        self.assert_call_counts(output, planned=2, attempted=3)
        self.assertEqual({
            row["requested_condition"]: float(row["raw_score"])
            for row in self.read_rows(output, "scores.tsv")
        }, {first: -0.125, second: 0.25})
