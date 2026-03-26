import os
import sys

import numpy as np


sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
)


def make_sample(bins_per_chr):
    sample = {}
    for i, n in enumerate(bins_per_chr, start=1):
        sample[str(i)] = np.full(int(n), 100.0 + i, dtype=float)
    return sample


def make_ref_arrays(bins_per_chr, refsize=3):
    n_bins = int(np.sum(bins_per_chr))
    mask = np.ones(n_bins, dtype=bool)
    masked_bins_per_chr = np.array([int(x) for x in bins_per_chr], dtype=int)
    masked_bins_per_chr_cum = np.cumsum(masked_bins_per_chr)

    n_components = 2
    pca_components = np.zeros((n_components, n_bins), dtype=float)
    pca_mean = np.ones(n_bins, dtype=float)

    indexes = np.zeros((n_bins, refsize), dtype=np.int32)
    distances = np.zeros((n_bins, refsize), dtype=float)

    max_valid = max(n_bins - 2, 0)
    base = np.arange(refsize, dtype=np.int32)
    base = np.clip(base, 0, max_valid)
    for i in range(n_bins):
        indexes[i, :] = base
        distances[i, :] = 0.1

    null_ratios = np.zeros((n_bins, min(10, n_bins)), dtype=float)

    return {
        "bins_per_chr": np.array(bins_per_chr, dtype=int),
        "mask": mask,
        "masked_bins_per_chr": masked_bins_per_chr,
        "masked_bins_per_chr_cum": masked_bins_per_chr_cum,
        "pca_components": pca_components,
        "pca_mean": pca_mean,
        "indexes": indexes,
        "distances": distances,
        "null_ratios": null_ratios,
    }
