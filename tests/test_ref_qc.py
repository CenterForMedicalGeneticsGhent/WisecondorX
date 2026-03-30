import unittest

import numpy as np

from wisecondorx.refqc import (
    _legacy_autosomal_prefix_compat_issues,
    _verdict_m,
    _get_gender_suffixes,
    _compute_per_bin_stats,
    _verdict_f,
)


class TestRefQc(unittest.TestCase):
    def test_get_gender_suffixes(self):
        self.assertEqual(_get_gender_suffixes({"bins_per_chr.F": 1}), [".F"])
        self.assertEqual(_get_gender_suffixes({"bins_per_chr.M": 1}), [".M"])
        self.assertEqual(_get_gender_suffixes({"bins_per_chr": 1}), [""])
        self.assertEqual(
            sorted(
                _get_gender_suffixes(
                    {"bins_per_chr.F": 1, "bins_per_chr.M": 1}
                )
            ),
            [".F", ".M"],
        )

    def test_compute_per_bin_stats_2d(self):
        indexes = np.zeros((2, 3))
        distances = np.array([[1, 2, 3], [4, 5, 6]])
        mean_d, max_d, n_refs = _compute_per_bin_stats(indexes, distances)
        np.testing.assert_array_equal(mean_d, [2, 5])
        np.testing.assert_array_equal(max_d, [3, 6])
        np.testing.assert_array_equal(n_refs, [3, 3])

    def test_compute_per_bin_stats_object(self):
        indexes = np.array([None, None], dtype=object)
        distances = np.array(
            [np.array([1, 2]), np.array([4, 5, 6])], dtype=object
        )
        mean_d, max_d, n_refs = _compute_per_bin_stats(indexes, distances)
        np.testing.assert_array_equal(mean_d, [1.5, 5])
        np.testing.assert_array_equal(max_d, [2, 6])
        np.testing.assert_array_equal(n_refs, [2, 3])

    def test_verdict_f_pass(self):
        m = {
            "n_valid": 10,
            "n_low_refs": 0,
            "std_of_means": 0.5,
            "outlier_pct": 0.5,
        }
        self.assertEqual(_verdict_f(m), ("PASS", ""))

    def test_verdict_f_fail_high_std(self):
        m = {
            "n_valid": 10,
            "n_low_refs": 0,
            "std_of_means": 15,
            "outlier_pct": 0.5,
        }
        v, msg = _verdict_f(m)
        self.assertEqual(v, "FAIL")
        self.assertIn("high", msg)

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
            "masked_bins_per_chr_cum.F": np.array(
                list(range(1, 23)) + [23], dtype=int
            ),
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
            "masked_bins_per_chr_cum.F": np.array(
                list(range(1, 23)) + [23], dtype=int
            ),
            "mask.F": np.array([True] * 23, dtype=bool),
        }
        issues = _legacy_autosomal_prefix_compat_issues(ref)
        self.assertEqual(issues, [])


if __name__ == "__main__":
    unittest.main()
