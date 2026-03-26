import unittest

import numpy as np

from wisecondorx.ref_qc import _legacy_autosomal_prefix_compat_issues, _verdict_m


class TestRefQc(unittest.TestCase):
    def test_ref_qc_male_verdict_warns_on_low_chry_usable_pct(self):
        metrics = {
            "n_valid": 1000,
            "n_low_refs": 0,
            "mean_of_means": 1.0,
            "outlier_pct": 0.0,
            "chrY": {
                "n_valid": 100,
                "mean_of_means": 1.0,
                "usable_pct": 40.0,
            },
        }
        verdict, msg = _verdict_m(metrics)
        self.assertEqual(verdict, "WARN")
        self.assertIn("usable chrY bins", msg)

    def test_ref_qc_male_verdict_fails_on_very_low_chry_usable_pct(self):
        metrics = {
            "n_valid": 1000,
            "n_low_refs": 0,
            "mean_of_means": 1.0,
            "outlier_pct": 0.0,
            "chrY": {
                "n_valid": 100,
                "mean_of_means": 1.0,
                "usable_pct": 10.0,
            },
        }
        verdict, msg = _verdict_m(metrics)
        self.assertEqual(verdict, "FAIL")
        self.assertIn("usable chrY bins", msg)

    def test_legacy_autosomal_prefix_compat_detects_mismatch(self):
        ref = {
            "bins_per_chr": np.array([1] * 22, dtype=int),
            "masked_bins_per_chr_cum": np.arange(1, 23, dtype=int),
            "mask": np.array([True] * 22, dtype=bool),
            "bins_per_chr.F": np.array([1] * 23, dtype=int),
            "masked_bins_per_chr_cum.F": np.array(list(range(1, 23)) + [23], dtype=int),
            "mask.F": np.array([True] * 21 + [False, True], dtype=bool),
        }
        issues = _legacy_autosomal_prefix_compat_issues(ref)
        self.assertTrue(len(issues) > 0)

    def test_legacy_autosomal_prefix_compat_accepts_matching_prefix(self):
        ref = {
            "bins_per_chr": np.array([1] * 22, dtype=int),
            "masked_bins_per_chr_cum": np.arange(1, 23, dtype=int),
            "mask": np.array([True] * 22, dtype=bool),
            "bins_per_chr.F": np.array([1] * 23, dtype=int),
            "masked_bins_per_chr_cum.F": np.array(list(range(1, 23)) + [23], dtype=int),
            "mask.F": np.array([True] * 23, dtype=bool),
        }
        issues = _legacy_autosomal_prefix_compat_issues(ref)
        self.assertEqual(issues, [])


if __name__ == "__main__":
    unittest.main()

