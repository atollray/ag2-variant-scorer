"""Offline checks for the multi-condition scoring CLI."""

import csv
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import grpc
import pandas as pd

import score


CONDITIONS = ("UBERON:0000053_CL:0000573", "UBERON:0000053_CL:0002586")
SCORERS = ("RNA_SEQ", "RNA_SEQ_ACTIVE")


def score_row(raw_score: float, imputation_mode: str) -> pd.DataFrame:
    return pd.DataFrame({
        "output_type": ["RNA_SEQ"],
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

    def run_cli(self, rows: str) -> int:
        self.variants.write_text(
            "input_id\tchromosome\tposition\tref\talt\n" + rows, encoding="utf-8"
        )
        args = [
            "score.py", "--input", str(self.variants), "--reference", str(self.reference),
            "--output-dir", str(self.output),
        ]
        for condition in CONDITIONS:
            args.extend(["--condition", condition])
        for scorer in SCORERS:
            args.extend(["--scorer", scorer])
        with mock.patch.object(sys, "argv", args), mock.patch.dict(
            os.environ, {"ALPHAGENOME_API_KEY": "offline-fixture-key"}
        ):
            return score.main()

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
