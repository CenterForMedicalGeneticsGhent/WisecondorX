import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from tests._test_utils import make_ref_arrays


class TestNewrefContract(unittest.TestCase):
    def test_tool_newref_merge_keeps_predict_contract_keys(self):
        from wisecondorx.newref import tool_newref_merge

        def _write_component(path, gender, n_chr):
            bins_per_chr = [1] * n_chr
            ref = make_ref_arrays(bins_per_chr, refsize=2)
            np.savez_compressed(
                path,
                binsize=100000,
                gender=gender,
                mask=ref["mask"],
                bins_per_chr=ref["bins_per_chr"],
                masked_bins_per_chr=ref["masked_bins_per_chr"],
                masked_bins_per_chr_cum=ref["masked_bins_per_chr_cum"],
                pca_components=ref["pca_components"],
                pca_mean=ref["pca_mean"],
                indexes=ref["indexes"],
                distances=ref["distances"],
                null_ratios=ref["null_ratios"],
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            a_path = os.path.join(tmpdir, "A.npz")
            f_path = os.path.join(tmpdir, "F.npz")
            m_path = os.path.join(tmpdir, "M.npz")
            out_path = os.path.join(tmpdir, "merged.npz")
            _write_component(a_path, "A", 22)
            _write_component(f_path, "F", 23)
            _write_component(m_path, "M", 24)

            args = SimpleNamespace(outfile=out_path, nipt=False)
            tool_newref_merge(
                args, [a_path, f_path, m_path], trained_cutoff=0.005
            )

            merged = np.load(out_path, allow_pickle=True)
            try:
                base_keys = [
                    "binsize",
                    "mask",
                    "bins_per_chr",
                    "masked_bins_per_chr",
                    "masked_bins_per_chr_cum",
                    "pca_components",
                    "pca_mean",
                    "indexes",
                    "distances",
                    "null_ratios",
                ]
                for k in base_keys:
                    self.assertIn(k, merged.keys())
                    self.assertIn("{}.F".format(k), merged.keys())
                    self.assertIn("{}.M".format(k), merged.keys())
                self.assertIn("has_female", merged.keys())
                self.assertIn("has_male", merged.keys())
                self.assertIn("is_nipt", merged.keys())
                self.assertIn("trained_cutoff", merged.keys())

                # Legacy predict compatibility invariants:
                # autosomal masked prefixes in .F/.M must match unsuffixed A.
                autosomal_masked_a = int(merged["masked_bins_per_chr_cum"][21])
                autosomal_span_a = int(np.sum(merged["bins_per_chr"][:22]))
                autosomal_prefix_a = merged["mask"][:autosomal_span_a]
                for suffix in [".F", ".M"]:
                    autosomal_masked_s = int(
                        merged["masked_bins_per_chr_cum{}".format(suffix)][21]
                    )
                    self.assertEqual(autosomal_masked_s, autosomal_masked_a)

                    autosomal_span_s = int(
                        np.sum(merged["bins_per_chr{}".format(suffix)][:22])
                    )
                    autosomal_prefix_s = merged["mask{}".format(suffix)][
                        :autosomal_span_s
                    ]
                    self.assertTrue(
                        np.array_equal(autosomal_prefix_s, autosomal_prefix_a)
                    )

                    # Old predict appends A autosomes + branch gonosomal tail and then
                    # inflates against branch mask; lengths must align.
                    total_masked_s = int(
                        merged["masked_bins_per_chr_cum{}".format(suffix)][-1]
                    )
                    gonosomal_tail_s = total_masked_s - autosomal_masked_s
                    appended_len = autosomal_masked_a + gonosomal_tail_s
                    self.assertEqual(
                        appended_len,
                        int(np.sum(merged["mask{}".format(suffix)])),
                    )
            finally:
                merged.close()

    def test_tool_newref_prep_male_pruning_keeps_autosomes_and_chry(self):
        from wisecondorx.newref import tool_newref_prep

        bins_per_chr = [10] * 24
        total_bins = int(np.sum(bins_per_chr))
        mask = np.ones(total_bins, dtype=bool)

        vals = np.resize(
            np.array([-3, -2, -1, 0, 1, 2, 3, -2, 2, 0], dtype=float),
            total_bins,
        )
        autosome_outlier = 5
        chry_start = int(np.sum(bins_per_chr[:23]))
        chry_outlier = chry_start + 5
        vals[autosome_outlier] = 7.0
        vals[chry_outlier] = 7.0
        pca_first = np.stack([vals, np.zeros(total_bins, dtype=float)], axis=1)

        call_count = {"n": 0}

        def fake_normalize_and_mask(samples, chrs, local_mask):
            return np.zeros((int(np.sum(local_mask)), 2), dtype=float)

        def fake_train_pca(ref_data):
            n_bins = ref_data.shape[0]
            pca_obj = SimpleNamespace(
                components_=np.zeros((1, n_bins), dtype=float),
                mean_=np.zeros(n_bins, dtype=float),
            )
            if call_count["n"] == 0:
                call_count["n"] += 1
                return pca_first, pca_obj
            return np.ones((n_bins, 2), dtype=float), pca_obj

        with tempfile.TemporaryDirectory() as tmpdir:
            args = SimpleNamespace(
                binsize=100000,
                prepdatafile=os.path.join(tmpdir, "prep_data.npy"),
                prepfile=os.path.join(tmpdir, "prep.npz"),
            )
            with (
                patch(
                    "wisecondorx.newref.normalize_and_mask",
                    side_effect=fake_normalize_and_mask,
                ),
                patch(
                    "wisecondorx.newref.train_pca",
                    side_effect=fake_train_pca,
                ),
            ):
                tool_newref_prep(
                    args,
                    samples=np.array([{}], dtype=object),
                    gender="M",
                    mask=mask,
                    bins_per_chr=bins_per_chr,
                )

            prep = np.load(args.prepfile, allow_pickle=True)
            try:
                saved_mask = prep["mask"]
                # In M branch, autosomes are inherited from A and must not be pruned further.
                self.assertTrue(saved_mask[autosome_outlier])
                self.assertTrue(saved_mask[chry_outlier])
            finally:
                prep.close()

    def test_tool_newref_prep_autosomal_pruning_still_applies_in_a_branch(
        self,
    ):
        from wisecondorx.newref import tool_newref_prep

        bins_per_chr = [10] * 22
        total_bins = int(np.sum(bins_per_chr))
        mask = np.ones(total_bins, dtype=bool)

        vals = np.resize(
            np.array([-3, -2, -1, 0, 1, 2, 3, -2, 2, 0], dtype=float),
            total_bins,
        )
        autosome_outlier = 5
        vals[autosome_outlier] = 7.0
        pca_first = np.stack([vals, np.zeros(total_bins, dtype=float)], axis=1)

        call_count = {"n": 0}

        def fake_normalize_and_mask(samples, chrs, local_mask):
            return np.zeros((int(np.sum(local_mask)), 2), dtype=float)

        def fake_train_pca(ref_data):
            n_bins = ref_data.shape[0]
            pca_obj = SimpleNamespace(
                components_=np.zeros((1, n_bins), dtype=float),
                mean_=np.zeros(n_bins, dtype=float),
            )
            if call_count["n"] == 0:
                call_count["n"] += 1
                return pca_first, pca_obj
            return np.ones((n_bins, 2), dtype=float), pca_obj

        with tempfile.TemporaryDirectory() as tmpdir:
            args = SimpleNamespace(
                binsize=100000,
                prepdatafile=os.path.join(tmpdir, "prep_data.npy"),
                prepfile=os.path.join(tmpdir, "prep.npz"),
            )
            with (
                patch(
                    "wisecondorx.newref.normalize_and_mask",
                    side_effect=fake_normalize_and_mask,
                ),
                patch(
                    "wisecondorx.newref.train_pca",
                    side_effect=fake_train_pca,
                ),
            ):
                tool_newref_prep(
                    args,
                    samples=np.array([{}], dtype=object),
                    gender="A",
                    mask=mask,
                    bins_per_chr=bins_per_chr,
                )

            prep = np.load(args.prepfile, allow_pickle=True)
            try:
                saved_mask = prep["mask"]
                self.assertFalse(saved_mask[autosome_outlier])
            finally:
                prep.close()


if __name__ == "__main__":
    unittest.main()
