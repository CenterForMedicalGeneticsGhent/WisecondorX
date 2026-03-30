import unittest
import os
import shutil
import tempfile
from pathlib import Path
import numpy as np
import argparse

from wisecondorx.newref import (
    _samples_to_matrix,
    get_mask,
    normalize_and_mask,
    train_gender_model,
    train_pca,
    get_ref_for_bins,
    _robust_cutoff,
    _masked_chr_ids,
    _split_by_chr,
    _get_part,
    force_remove,
    tool_newref_merge,
    tool_newref_prep,
)


class TestNewref(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        # Create dummy samples
        # Each sample is a dict of arrays for chromosomes '1' to '24'
        self.bins_per_chr = [10] * 24
        self.samples = []
        for i in range(10):
            sample = {}
            for c in range(1, 25):
                # Add some variation
                sample[str(c)] = np.random.poisson(
                    100, self.bins_per_chr[c - 1]
                ).astype(float)
            self.samples.append(sample)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_samples_to_matrix(self):
        matrix = _samples_to_matrix(self.samples, ["1", "2"])
        self.assertEqual(matrix.shape, (20, 10))  # (10+10 bins, 10 samples)
        self.assertTrue(np.all(matrix >= 0))

    def test_get_mask(self):
        # Introduce a zero bin
        self.samples[0]["1"][0] = 0
        self.samples[1]["1"][0] = 0
        mask, bins = get_mask(self.samples)
        self.assertEqual(len(mask), sum(self.bins_per_chr))
        self.assertEqual(len(bins), 24)
        # The first bin should likely be masked if it's 0 in most or low
        # In our case, 2/10 are zero, median might still be high enough.
        # Let's force it to be masked by making it zero everywhere.
        for s in self.samples:
            s["1"][1] = 0
        mask, bins = get_mask(self.samples)
        self.assertFalse(mask[1])

    def test_normalize_and_mask(self):
        mask = np.ones(sum(self.bins_per_chr), dtype=bool)
        mask[0] = False
        masked_data = normalize_and_mask(
            self.samples, [str(i) for i in range(1, 25)], mask
        )
        self.assertEqual(masked_data.shape[0], sum(self.bins_per_chr) - 1)
        self.assertEqual(masked_data.shape[1], 10)

    def test_train_gender_model(self):
        # Create samples with distinct Y fractions
        gender_samples = []
        for i in range(10):
            s = {str(c): np.ones(10) * 100 for c in range(1, 25)}
            if i < 5:  # "Females" - low Y
                s["24"] = np.ones(10) * 1
            else:  # "Males" - high Y
                s["24"] = np.ones(10) * 10
            gender_samples.append(s)

        args = argparse.Namespace(plotyfrac=None, yfrac=None)
        genders, cutoff = train_gender_model(args, gender_samples)
        self.assertEqual(genders.count("F"), 5)
        self.assertEqual(genders.count("M"), 5)
        self.assertTrue(0 < cutoff < 0.01)

    def test_train_pca(self):
        data = np.random.rand(100, 20)  # 100 bins, 20 samples
        corrected, pca = train_pca(data, pcacomp=2)
        self.assertEqual(corrected.shape, (100, 20))
        self.assertEqual(len(pca.components_), 2)

    def test_get_ref_for_bins(self):
        ref_size = 5
        pca_corrected_data = np.random.rand(50, 10)
        chr_data = np.random.rand(100, 10)
        ref_idx, ref_dist = get_ref_for_bins(
            ref_size, 0, 10, pca_corrected_data, chr_data
        )
        self.assertEqual(ref_idx.shape, (10, ref_size))
        self.assertEqual(ref_dist.shape, (10, ref_size))
        # Distances should be sorted
        for i in range(10):
            self.assertTrue(np.all(np.diff(ref_dist[i]) >= 0))

    def test_robust_cutoff(self):
        distances = np.array([1, 1.1, 1.2, 1.3, 1.4, 10.0])  # 10 is outlier
        cutoff = _robust_cutoff(distances, mad_multiplier=3, floor=0)
        self.assertTrue(cutoff > 1.4)
        self.assertTrue(cutoff < 10.0)

    def test_masked_chr_ids(self):
        mask = np.array([True, False, True, True])
        bins_per_chr = [2, 2]
        indices, ids = _masked_chr_ids(mask, bins_per_chr)
        # indices should be [0, 2, 3]
        # ids should be [1, 2, 2]
        np.testing.assert_array_equal(indices, [0, 2, 3])
        np.testing.assert_array_equal(ids, [1, 2, 2])

    def test_split_by_chr(self):
        chr_bin_sums = [10, 25, 40]
        # start=5, end=30
        areas = _split_by_chr(5, 30, chr_bin_sums)
        # Expected:
        # Area 1: chr 0, start 5, end 10
        # Area 2: chr 1, start 10, end 25
        # Area 3: chr 2, start 25, end 30
        self.assertEqual(len(areas), 3)
        self.assertEqual(areas[0], [0, 5, 10])
        self.assertEqual(areas[1], [1, 10, 25])
        self.assertEqual(areas[2], [2, 25, 30])

    def test_get_part(self):
        start, end = _get_part(0, 2, 100)
        self.assertEqual((start, end), (0, 50))
        start, end = _get_part(1, 2, 100)
        self.assertEqual((start, end), (50, 100))

    def test_force_remove(self):
        tmp_file = os.path.join(self.test_dir, "temp.txt")
        Path(tmp_file).touch()
        self.assertTrue(os.path.exists(tmp_file))
        force_remove(tmp_file)
        self.assertFalse(os.path.exists(tmp_file))

    def test_tool_newref_merge(self):
        # Create dummy partial npz files
        outfiles = []
        for g in ["A", "F"]:
            p = os.path.join(self.test_dir, f"part_{g}.npz")
            outfiles.append(p)
            np.savez_compressed(p, gender=g, data=np.array([1, 2, 3]))

        outfile = os.path.join(self.test_dir, "final.npz")
        args = argparse.Namespace(outfile=outfile, nipt=False)
        tool_newref_merge(args, outfiles, trained_cutoff=0.005)

        self.assertTrue(os.path.exists(outfile))
        final = np.load(outfile, allow_pickle=True)
        self.assertTrue("data" in final)  # From A
        self.assertTrue("data.F" in final)
        self.assertEqual(final["trained_cutoff"], 0.005)

    def test_tool_newref_prep(self):
        # This function writes files, so we test if it creates them and their content
        args = argparse.Namespace(
            binsize=5000,
            prepfile=os.path.join(self.test_dir, "prep.npz"),
            prepdatafile=os.path.join(self.test_dir, "prep_data.npy"),
        )
        mask = np.ones(sum(self.bins_per_chr), dtype=bool)
        tool_newref_prep(args, self.samples, "A", mask, self.bins_per_chr)

        self.assertTrue(os.path.exists(args.prepfile))
        self.assertTrue(os.path.exists(args.prepdatafile))

        prep = np.load(args.prepfile, allow_pickle=True)
        self.assertEqual(prep["gender"], "A")
        self.assertEqual(len(prep["bins_per_chr"]), 22)


if __name__ == "__main__":
    unittest.main()
