"""Offline checks for the multi-condition scoring CLI."""

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

import grpc
import pandas as pd

from ag2_variant_scorer import score


CONDITIONS = ("UBERON:0000053_CL:0000573", "UBERON:0000053_CL:0002586")
SCORERS = ("RNA_SEQ", "RNA_SEQ_ACTIVE")


def score_row(raw_score: float, imputation_mode: str, output_type: str = "RNA_SEQ") -> pd.DataFrame:
    return pd.DataFrame({
        "output_type": [output_type],
        "raw_score": [raw_score],
        "imputation_mode": [imputation_mode],
    })


class Unavailable(grpc.RpcError):
    def code(self):
        return grpc.StatusCode.UNAVAILABLE

    def details(self):
        return "Offline fixture RPC failure."


class ScoreCliTest(unittest.TestCase):
    def setUp(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        self.reference = directory / "reference.fa"
        self.variants = directory / "variants.tsv"
        self.output = directory / "output"
        bases = b"A" * 1_200_000
        self.reference.write_bytes(b">chr22\n" + bases + b"\n")
        self.reference.with_suffix(".fa.fai").write_bytes(
            f"chr22\t{len(bases)}\t7\t{len(bases)}\t{len(bases) + 1}\n".encode("ascii")
        )
        self.model = self.enterContext(mock.patch.object(score.dna_client, "create")).return_value
        self.model.supported_conditions.return_value = pd.DataFrame({
            "condition": CONDITIONS,
            "biosample_name": ["macula cone cell", "macula retinal pigment epithelial cell"],
            "rna_seq": ["imputed_supported", "imputed_supported"],
        })
        self.enterContext(mock.patch.object(
            score.variant_scorers, "tidy_anndata", side_effect=lambda table, **_: table.copy()
        ))

    def run_cli(self, rows: str, scorers: tuple[str, ...] = SCORERS) -> int:
        options = ["--output-dir", str(self.output)]
        for condition in CONDITIONS:
            options.extend(["--condition", condition])
        for scorer in scorers:
            options.extend(["--scorer", scorer])
        return self.invoke_cli(rows, options)

    def invoke_cli(self, rows: str, options: list[str]) -> int:
        self.variants.write_text(
            "input_id\tchromosome\tposition\tref\talt\n" + rows, encoding="utf-8"
        )
        args = [
            "score.py", "--input", str(self.variants), "--reference", str(self.reference),
            *options,
        ]
        with mock.patch.object(sys, "argv", args), mock.patch.dict(
            os.environ, {"ALPHAGENOME_API_KEY": "offline-fixture-key"}
        ):
            return score.main()

    def write_panel(self, condition_scorers: dict[str, list[str]]) -> Path:
        path = self.reference.with_name("panel.toml")
        condition_lines = [
            f"{json.dumps(condition)} = {json.dumps(names)}"
            for condition, names in condition_scorers.items()
        ]
        lines = ['name = "test panel"', 'version = 1', '[conditions]', *condition_lines]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def read_tsv(self, name: str) -> list[dict[str, str]]:
        with (self.output / name).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_group_scorers_per_condition_and_retain_their_names(self):
        self.model.score_variant.side_effect = [
            [score_row(-0.5, "imputed_supported"), score_row(8.0, "imputed_supported")],
            [score_row(0.25, "imputed_supported"), score_row(3.0, "imputed_supported")],
        ]

        self.assertEqual(
            self.run_cli("checked\tchr22\t600001\tA\tC\nwrong_ref\tchr22\t600002\tC\tA\n"), 1
        )

        calls = [call.kwargs for call in self.model.score_variant.call_args_list]
        self.assertEqual([call["condition"] for call in calls], list(CONDITIONS))
        for call in calls:
            self.assertEqual(call["variant"].position, 600001)
            self.assertEqual(list(call["variant_scorers"]), [
                score.variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in SCORERS
            ])
        self.model.predict_variant.assert_not_called()
        self.assertFalse((self.output / "splicing_summary.tsv").exists())

        rows = self.read_tsv("scores.tsv")
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["input_id"] for row in rows}, {"checked"})
        self.assertEqual({row["output_type"] for row in rows}, {"RNA_SEQ"})
        self.assertEqual({row["imputation_mode"] for row in rows}, {"imputed_supported"})
        self.assertEqual({
            (row["requested_condition"], row["requested_scorer"]): float(row["raw_score"])
            for row in rows
        }, {
            (CONDITIONS[0], SCORERS[0]): -0.5,
            (CONDITIONS[0], SCORERS[1]): 8.0,
            (CONDITIONS[1], SCORERS[0]): 0.25,
            (CONDITIONS[1], SCORERS[1]): 3.0,
        })
        self.assertEqual([row["input_id"] for row in self.read_tsv("rejected.tsv")], ["wrong_ref"])

    def test_empty_and_failed_requests_are_visible_in_the_ledger(self):
        self.model.score_variant.side_effect = [
            [pd.DataFrame(), score_row(2.0, "imputed_supported")],
            Unavailable("Offline fixture RPC failure."),
        ]

        self.assertEqual(self.run_cli("checked\tchr22\t600001\tA\tC\n"), 1)

        scores = self.read_tsv("scores.tsv")
        self.assertEqual(len(scores), 1)
        self.assertEqual(scores[0]["requested_scorer"], SCORERS[1])
        requests = self.read_tsv("requests.tsv")
        by_pair = {
            (row["requested_condition"], row["requested_scorer"]): (
                row["status"], int(row["score_rows"])
            )
            for row in requests
        }
        self.assertEqual(by_pair, {
            (CONDITIONS[0], SCORERS[0]): ("no_rows", 0),
            (CONDITIONS[0], SCORERS[1]): ("scored", 1),
            (CONDITIONS[1], SCORERS[0]): ("api_error", 0),
            (CONDITIONS[1], SCORERS[1]): ("api_error", 0),
        })

    def test_splicing_summary_uses_observed_components_per_condition(self):
        sites, usage, junctions = ("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
        self.model.score_variant.side_effect = [
            [
                pd.concat([
                    score_row(-0.25, "observed", sites),
                    score_row(0.95, "imputed_exploratory", sites),
                ], ignore_index=True),
                score_row(0.5, "observed", usage),
                score_row(4.0, "observed", junctions),
            ],
            [
                score_row(0.4, "observed", sites),
                score_row(0.7, "imputed_exploratory", usage),
                score_row(2.0, "observed", junctions),
            ],
            Unavailable("Offline fixture RPC failure."),
        ]

        self.assertEqual(self.run_cli(
            "checked\tchr22\t600001\tA\tC\ninterrupted\tchr22\t600002\tA\tT\n",
            scorers=(sites, usage, junctions),
        ), 1)

        rows = self.read_tsv("splicing_summary.tsv")
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["input_id"] for row in rows}, {"checked"})
        by_condition = {row["requested_condition"]: row for row in rows}
        full = by_condition[CONDITIONS[0]]
        self.assertEqual([float(full[name]) for name in (sites, usage, junctions)], [0.25, 0.5, 4.0])
        self.assertAlmostEqual(float(full["merged_splicing"]), 1.55)
        self.assertEqual(full["missing_observed_components"], "")

        incomplete = by_condition[CONDITIONS[1]]
        self.assertEqual(incomplete[usage], "")
        self.assertEqual(incomplete["merged_splicing"], "")
        self.assertEqual(incomplete["missing_observed_components"], usage)

    def test_panel_dry_run_checks_ref_without_model_or_output(self):
        panel = Path(__file__).with_name("panels") / "retina-rpe.toml"
        stdout = io.StringIO()
        with redirect_stdout(stdout), mock.patch.object(score, "api_key_from_file_or_environment") as key:
            status = self.invoke_cli(
                "checked\tchr22\t600001\tA\tC\nwrong_ref\tchr22\t600002\tC\tA\n",
                ["--panel", str(panel), "--dry-run"],
            )
        self.assertEqual(status, 1)
        key.assert_not_called()
        score.dna_client.create.assert_not_called()
        self.assertFalse(self.output.exists())
        self.assertEqual(stdout.getvalue().splitlines(), [
            "REF-verified rows: 1", "Rejected rows: 1", "Conditions: 4", "Planned score calls: 4",
        ])

    def test_panel_records_uneven_scorer_sets_and_preserves_results_after_error(self):
        first, second = CONDITIONS
        sites, usage, junctions = score.SPLICE_OUTPUTS
        condition_scorers = {
            first: ["SC_GENE_EXPRESSION"],
            second: [sites, usage, junctions],
        }
        panel = self.write_panel(condition_scorers)
        self.model.score_variant.side_effect = [
            [score_row(0.2, "observed", "SC_GENE_EXPRESSION")],
            [score_row(0.25, "observed", name) for name in condition_scorers[second]],
            Unavailable("Offline fixture RPC failure."),
        ]
        self.assertEqual(self.invoke_cli(
            "checked\tchr22\t600001\tA\tC\ninterrupted\tchr22\t600002\tA\tT\n",
            ["--panel", str(panel), "--output-dir", str(self.output)],
        ), 1)

        calls = [call.kwargs for call in self.model.score_variant.call_args_list]
        self.assertEqual([call["condition"] for call in calls], [first, second, first])
        for call in calls:
            self.assertEqual(list(call["variant_scorers"]), [
                score.variant_scorers.RECOMMENDED_VARIANT_SCORERS[name]
                for name in condition_scorers[call["condition"]]
            ])
        manifest = json.loads((self.output / "request.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["planned_score_calls"], 4)
        self.assertEqual(manifest["attempted_score_calls"], 3)
        self.assertEqual(manifest["condition_scorers"], condition_scorers)
        self.assertEqual(manifest["panel"], {
            "name": "test panel", "version": 1, "conditions": condition_scorers,
        })
        self.assertEqual(manifest["splicing_summary_tracks"], "observed")
        self.assertEqual([row["requested_condition"] for row in self.read_tsv("splicing_summary.tsv")], [second])
        self.assertEqual([row["status"] for row in self.read_tsv("requests.tsv")], ["scored"] * 4 + ["api_error"])
