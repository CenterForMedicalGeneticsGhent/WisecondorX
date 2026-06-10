import unittest
import numpy as np
import pandas as pd
from wisecondorx.predict_tools import get_segment_zscore


class TestGetSegmentZScore(unittest.TestCase):
    def test_basic_zscore_calculation(self):
        # Create a segments dataframe
        segments_df = pd.DataFrame(
            {
                "sample": ["test_sample"],
                "chrom": ["1"],
                "loc.start": [1],
                "loc.end": [3],
                "num_mark": [3],
                "seg.mean": [1.2],
            }
        )

        # Create cna_df with 3 bins
        cna_df = pd.DataFrame(
            {
                "chrom": ["1", "1", "1"],
                "maploc": [1, 2, 3],
                "null_ratios_mean": [1.0, 1.0, 1.0],
                "null_ratios_sd": [0.1, 0.1, 0.1],
                "ratios": [1.2, 1.3, 1.1],
                "weights": [1.0, 1.0, 1.0],
                "zscores": [2.0, 3.0, 1.0],
            }
        )

        res_df = get_segment_zscore(segments_df, cna_df)

        # Expected null_mean = 1.0
        # Expected null_sd = sqrt(3 * 0.1^2) / 3 = sqrt(0.03) / 3 = 0.05773502691896258
        # Expected zscore = (1.2 - 1.0) / 0.05773502691896258 = 3.4641016151377544
        self.assertIn("zscore", res_df.columns)
        self.assertAlmostEqual(
            res_df.loc[0, "zscore"], 3.4641016151377544, places=5
        )

    def test_filtering_zero_ratios(self):
        # One bin has ratio = 0, which should be filtered out
        segments_df = pd.DataFrame(
            {
                "sample": ["test_sample"],
                "chrom": ["1"],
                "loc.start": [1],
                "loc.end": [3],
                "num_mark": [3],
                "seg.mean": [1.2],
            }
        )

        # Bin at maploc 2 has ratios = 0, meaning it should be filtered out
        cna_df = pd.DataFrame(
            {
                "chrom": ["1", "1", "1"],
                "maploc": [1, 2, 3],
                "null_ratios_mean": [1.0, 5.0, 1.0],  # 5.0 should be ignored
                "null_ratios_sd": [0.1, 0.5, 0.1],  # 0.5 should be ignored
                "ratios": [1.2, 0.0, 1.1],  # 0.0 is bad, should filter out
                "weights": [1.0, 1.0, 1.0],
                "zscores": [2.0, 0.0, 1.0],
            }
        )

        res_df = get_segment_zscore(segments_df, cna_df)

        # Only bin 1 and 3 are valid
        # Expected null_mean = (1.0 * 1.0 + 1.0 * 1.0) / 2.0 = 1.0
        # Expected null_sd = sqrt(0.1^2 * 1.0^2 + 0.1^2 * 1.0^2) / 2.0 = sqrt(0.02) / 2.0 = 0.07071067811865475
        # Expected zscore = (1.2 - 1.0) / 0.07071067811865475 = 2.8284271247461903
        self.assertAlmostEqual(
            res_df.loc[0, "zscore"], 2.8284271247461903, places=5
        )

    def test_infinite_values_cleaned(self):
        segments_df = pd.DataFrame(
            {
                "sample": ["test_sample"],
                "chrom": ["1"],
                "loc.start": [1],
                "loc.end": [3],
                "num_mark": [3],
                "seg.mean": [1.2],
            }
        )

        # One bin has infinite null_ratios_mean, which should be replaced by NaN and thus ignored
        cna_df = pd.DataFrame(
            {
                "chrom": ["1", "1", "1"],
                "maploc": [1, 2, 3],
                "null_ratios_mean": [1.0, np.inf, 1.0],
                "null_ratios_sd": [0.1, 0.1, 0.1],
                "ratios": [1.2, 1.3, 1.1],
                "weights": [1.0, 1.0, 1.0],
                "zscores": [2.0, 3.0, 1.0],
            }
        )

        res_df = get_segment_zscore(segments_df, cna_df)

        # Expected null_mean = 1.0 (from bins 1 and 3)
        # Expected null_sd = sqrt(0.02) / 2.0 = 0.07071067811865475
        # Expected zscore = (1.2 - 1.0) / 0.07071067811865475 = 2.8284271247461903
        self.assertAlmostEqual(
            res_df.loc[0, "zscore"], 2.8284271247461903, places=5
        )

    def test_zero_null_sd_returns_nan(self):
        segments_df = pd.DataFrame(
            {
                "sample": ["test_sample"],
                "chrom": ["1"],
                "loc.start": [1],
                "loc.end": [3],
                "num_mark": [3],
                "seg.mean": [1.2],
            }
        )

        # null_ratios_sd is 0, so null_sd will be 0. Should return "nan"
        cna_df = pd.DataFrame(
            {
                "chrom": ["1", "1", "1"],
                "maploc": [1, 2, 3],
                "null_ratios_mean": [1.0, 1.0, 1.0],
                "null_ratios_sd": [0.0, 0.0, 0.0],
                "ratios": [1.2, 1.3, 1.1],
                "weights": [1.0, 1.0, 1.0],
                "zscores": [2.0, 3.0, 1.0],
            }
        )

        res_df = get_segment_zscore(segments_df, cna_df)
        self.assertEqual(res_df.loc[0, "zscore"], "nan")


if __name__ == "__main__":
    unittest.main()
