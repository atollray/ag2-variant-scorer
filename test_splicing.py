import unittest

import pandas as pd
from alphagenome.data.genome import Interval

from ag2_variant_scorer.splicing import junction_row, observed_summary


class SplicingTest(unittest.TestCase):
    def test_junction_intervals_are_oriented_by_strand(self):
        self.assertEqual(junction_row(Interval("chr22", 10, 20, strand="+"))[-2:], (10, 20))
        self.assertEqual(junction_row(Interval("chr22", 10, 20, strand="-"))[-2:], (20, 10))
        self.assertEqual(junction_row(Interval("chr22", 10, 20, strand="."))[-2:], (None, None))

    def test_merged_score_uses_only_finite_observed_components(self):
        rows = pd.DataFrame(
            [
                ("SPLICE_SITES", "observed", 0.25),
                ("SPLICE_SITE_USAGE", "observed", 0.5),
                ("SPLICE_JUNCTIONS", "observed", 4.0),
                ("SPLICE_JUNCTIONS", "observed", float("inf")),
                ("SPLICE_JUNCTIONS", "imputed_supported", 40.0),
            ],
            columns=["output_type", "imputation_mode", "raw_score"],
        )
        self.assertAlmostEqual(observed_summary(rows)["merged_splicing"], 1.55)

        incomplete = observed_summary(rows.drop(index=2))
        self.assertIsNone(incomplete["merged_splicing"])
        self.assertEqual(incomplete["missing_observed_components"], "SPLICE_JUNCTIONS")
