import logging
import sys
from enum import Enum
from typing import Dict

import numpy as np


class Sex(Enum):
    MALE = "M"
    FEMALE = "F"


def scale_sample(
    sample: Dict[str, np.ndarray], source_binsize: float, target_binsize: float
) -> Dict[str, np.ndarray]:
    """
    Scales the bin size of a sample.npz to the one
    requested for the reference
    """
    if source_binsize == target_binsize:
        return sample

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

    scaled_sample: Dict[str, np.ndarray] = dict()
    scale = target_binsize / source_binsize
    for chr_name in sample:
        chr_data = sample[chr_name]
        new_len = int(np.ceil(len(chr_data) / float(scale)))
        scaled_chr = np.zeros(new_len, dtype=np.int32)
        for i in range(new_len):
            scaled_chr[i] = np.sum(
                chr_data[int(i * scale) : int(i * scale + scale)]
            )
            scaled_sample[chr_name] = scaled_chr
    return scaled_sample


def sex_correct(
    sample: Dict[str, np.ndarray], sex: Sex
) -> Dict[str, np.ndarray]:
    """
    Levels gonosomal reads with the one at the autosomes.
    """
    if sex == Sex.MALE:
        sample["23"] = sample["23"] * 2
        sample["24"] = sample["24"] * 2

    return sample
