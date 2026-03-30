import unittest
import numpy as np
from wisecondorx.predict import (
    predict_sex,
    coverage_normalize_and_mask,
    project_pc,
    get_optimal_cutoff,
    get_weights,
    inflate_results,
    log_trans,
    _get_aberration_cutoff,
)


class TestPredict(unittest.TestCase):
    def test_predict_sex(self):
        sample = {"24": np.array([100]), "1": np.array([1000])}
        # Y fraction = 100 / 1100 = 0.09
        self.assertEqual(predict_sex(sample, 0.05), "M")
        self.assertEqual(predict_sex(sample, 0.15), "F")

    def test_coverage_normalize_and_mask(self):
        sample = {"1": np.array([10, 20]), "2": np.array([30])}
        ref_file = {
            "bins_per_chr": np.array([2, 1]),
            "mask": np.array([True, False, True]),
        }
        # Total reads = 10+20+30 = 60
        # Norm: [10/60, 20/60, 30/60] = [1/6, 1/3, 1/2]
        # Masked: [1/6, 1/2]
        masked = coverage_normalize_and_mask(sample, ref_file, "")
        np.testing.assert_allclose(masked, [1 / 6, 1 / 2])

    def test_project_pc(self):
        sample_data = np.array([1.0, 2.0, 3.0])
        pca_components = np.array([[0.1, 0.2, 0.3]])
        pca_mean = np.array([0.5, 0.5, 0.5])
        ref_file = {"pca_components": pca_components, "pca_mean": pca_mean}
        # centered = [0.5, 1.5, 2.5]
        # transform = 0.5*0.1 + 1.5*0.2 + 2.5*0.3 = 0.05 + 0.3 + 0.75 = 1.1
        # reconstructed = 1.1 * [0.1, 0.2, 0.3] + [0.5, 0.5, 0.5] = [0.11+0.5, 0.22+0.5, 0.33+0.5] = [0.61, 0.72, 0.83]
        # corrected = [1/0.61, 2/0.72, 3/0.83]
        corrected = project_pc(sample_data, ref_file, "")
        expected = sample_data / (
            np.dot(
                np.dot(sample_data - pca_mean, pca_components.T),
                pca_components,
            )
            + pca_mean
        )
        np.testing.assert_allclose(corrected, expected)

    def test_get_optimal_cutoff(self):
        ref_file = {"distances": np.array([1, 1, 1, 1, 100])}
        # average of [1,1,1,1,100] = 20.8
        # stddev is large.
        # Repeating 1 cycle:
        # average = 20.8, std = ~40, cutoff = 20.8 + 120 = 140.8
        # Cycle 2: same.
        # If we use a smaller set: [1, 1.1, 1.2, 10]
        ref_file = {"distances": np.array([1, 1.1, 1.2, 10])}
        cutoff = get_optimal_cutoff(ref_file, 2)
        self.assertTrue(cutoff > 1.2)

    def test_get_weights(self):
        ref_file = {"distances": np.array([4, 9, 16])}
        # sqrt: [2, 3, 4]
        # inverse_weights: [2, 3, 4] (because np.mean of a single value is the value)
        # weights: [1/2, 1/3, 1/4]
        # normalized: [0.5, 0.33, 0.25] / mean(0.5, 0.33, 0.25)
        # Wait, get_weights in code:
        # inverse_weights = [np.mean(np.sqrt(x)) for x in ref_file["distances"]]
        # This assumes ref_file["distances"] is a list of arrays?
        # Looking at newref.py: distances is an array of shape (bins, ref_size)
        distances = np.array([[4, 4], [9, 9]])
        ref_file = {"distances": distances}
        # sqrt: [[2, 2], [3, 3]]
        # inverse_weights: [2, 3]
        # weights: [1/2, 1/3] = [0.5, 0.333]
        # mean: 0.41666
        # normalized: [1.2, 0.8]
        weights = get_weights(ref_file, "")
        np.testing.assert_allclose(
            weights,
            [0.5 / np.mean([0.5, 1 / 3]), (1 / 3) / np.mean([0.5, 1 / 3])],
        )

    def test_inflate_results(self):
        results = np.array([10, 20])
        rem_input = {"mask": [True, False, True, False]}
        inflated = inflate_results(results, rem_input)
        self.assertEqual(inflated, [10, 0, 20, 0])

    def test_log_trans(self):
        results = {
            "results_r": [np.array([1.0, 2.0]), np.array([4.0])],
            "results_z": [np.array([0.1, 0.2]), np.array([0.3])],
            "results_w": [np.array([1.0, 1.0]), np.array([1.0])],
        }
        # log2: [[0, 1], [2]]
        # log_r_median = 1
        # final: [[-1, 0], [1]]
        log_trans(results, 1.0)
        np.testing.assert_allclose(results["results_r"][0], [-1, 0])
        np.testing.assert_allclose(results["results_r"][1], [1])

    def test_get_aberration_cutoff(self):
        # ploidy 2, beta 0.5
        # loss: log2((2 - 0.25)/2) = log2(1.75/2) = log2(0.875)
        # gain: log2((2 + 0.25)/2) = log2(2.25/2) = log2(1.125)
        loss, gain = _get_aberration_cutoff(0.5, 2)
        self.assertAlmostEqual(loss, np.log2(0.875))
        self.assertAlmostEqual(gain, np.log2(1.125))


if __name__ == "__main__":
    unittest.main()
