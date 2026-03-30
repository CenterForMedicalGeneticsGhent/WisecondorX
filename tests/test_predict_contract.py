import unittest

import numpy as np

from tests._test_utils import make_ref_arrays, make_sample


class TestPredictContract(unittest.TestCase):
    def test_predict_normalize_accepts_reference_contract_autosomes(self):
        from wisecondorx.predict import normalize

        bins_per_chr_a = [1] * 22
        ref = make_ref_arrays(bins_per_chr_a, refsize=3)
        sample = make_sample([1] * 24)

        class Args:
            maskrepeats = 1

        results_r, results_z, results_w, ref_sizes, m_lr, m_z = normalize(
            Args(), sample, ref, "A"
        )

        expected = int(np.sum(bins_per_chr_a))
        self.assertEqual(len(results_r), expected)
        self.assertEqual(len(results_z), expected)
        self.assertEqual(len(results_w), expected)
        self.assertEqual(len(ref_sizes), expected)
        self.assertEqual(np.ndim(results_r), 1)
        self.assertEqual(np.ndim(results_z), 1)

    def test_predict_normalize_accepts_reference_contract_male_gonosomes(self):
        from wisecondorx.predict import normalize

        bins_per_chr_m = [1] * 24
        ref_m = make_ref_arrays(bins_per_chr_m, refsize=3)
        ref = {
            # get_optimal_cutoff() reads unsuffixed distances as global cutoff input
            "distances": ref_m["distances"],
            "bins_per_chr.M": ref_m["bins_per_chr"],
            "mask.M": ref_m["mask"],
            "masked_bins_per_chr.M": ref_m["masked_bins_per_chr"],
            "masked_bins_per_chr_cum.M": ref_m["masked_bins_per_chr_cum"],
            "pca_components.M": ref_m["pca_components"],
            "pca_mean.M": ref_m["pca_mean"],
            "indexes.M": ref_m["indexes"],
            "distances.M": ref_m["distances"],
            "null_ratios.M": ref_m["null_ratios"],
        }
        sample = make_sample([1] * 24)

        class Args:
            maskrepeats = 1

        results_r, results_z, results_w, ref_sizes, m_lr, m_z = normalize(
            Args(), sample, ref, "M"
        )

        self.assertEqual(len(results_r), 2)
        self.assertEqual(len(results_z), 2)
        self.assertEqual(len(results_w), 2)
        self.assertEqual(len(ref_sizes), 2)
        self.assertEqual(np.ndim(results_r), 1)
        self.assertEqual(np.ndim(results_z), 1)

    def test_legacy_append_inflate_contract_requires_aligned_lengths(self):
        from wisecondorx.predict import inflate_results

        # Old predict behavior: len(appended_results) must equal number of True bins in
        # branch mask passed to inflate_results.
        rem_input = {"mask": np.array([True, False, True, True, False, True])}
        appended_results = np.array([1.0, 2.0, 3.0, 4.0])

        inflated = inflate_results(appended_results, rem_input)
        self.assertEqual(len(inflated), len(rem_input["mask"]))
        self.assertEqual(
            sum(x != 0 for x in inflated), int(np.sum(rem_input["mask"]))
        )


if __name__ == "__main__":
    unittest.main()
