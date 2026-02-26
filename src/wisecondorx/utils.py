import logging
import sys
from enum import Enum

import numpy as np


class Sex(Enum):
    MALE = "M"
    FEMALE = "F"
    AUTOSOMAL = "A"


def scale_bins_per_chromosome(
    bins_per_chr: dict[str, np.ndarray],
    source_binsize: float,
    target_binsize: float,
) -> dict[str, np.ndarray]:
    """
    Scales the bin size of a sample bins per chromosome to the one requested for the reference
    """
    if source_binsize == target_binsize:
        return bins_per_chr

    if (
        target_binsize == 0
        or source_binsize == 0
        or target_binsize < source_binsize
        or target_binsize % source_binsize > 0
    ):
        logging.critical(
            f"Impossible binsize scaling requested: {source_binsize} to {target_binsize}."
        )
        sys.exit(1)

    scaled_bins_per_chr: dict[str, np.ndarray] = dict()
    scale = target_binsize / source_binsize
    for chr_name in bins_per_chr:
        chr_data = bins_per_chr[chr_name]
        new_len = int(np.ceil(len(chr_data) / float(scale)))
        scaled_chr = np.zeros(new_len, dtype=np.int32)
        for i in range(new_len):
            scaled_chr[i] = np.sum(
                chr_data[int(i * scale) : int(i * scale + scale)]
            )
            scaled_bins_per_chr[chr_name] = scaled_chr
    return scaled_bins_per_chr


def sex_correct(
    sample: dict[str, np.ndarray], sex: Sex
) -> dict[str, np.ndarray]:
    """
    Levels gonosomal reads with the one at the autosomes.
    """
    if sex == Sex.MALE:
        sample["23"] = sample["23"] * 2
        sample["24"] = sample["24"] * 2

    return sample
