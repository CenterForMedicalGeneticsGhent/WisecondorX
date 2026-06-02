# WisecondorX

import logging
import sys
import math

import numpy as np
import pandas as pd
from enum import Enum


class Sex(Enum):
    MALE = "M"
    FEMALE = "F"
    AUTOSOMAL = "A"


def scale_sample(sample, from_size, to_size):
    """
    Scales the bin size of a sample.npz to the one
    requested for the reference
    """
    if not to_size or from_size == to_size:
        return sample

    if (
        to_size == 0
        or from_size == 0
        or to_size < from_size
        or to_size % from_size > 0
    ):
        logging.critical(
            "Impossible binsize scaling requested: {} to {}".format(
                int(from_size), int(to_size)
            )
        )
        sys.exit()

    return_sample = dict()
    scale = to_size / from_size
    for chr_name in sample:
        chr_data = sample[chr_name]
        new_len = int(np.ceil(len(chr_data) / float(scale)))
        scaled_chr = np.zeros(new_len, dtype=np.int32)
        for i in range(new_len):
            scaled_chr[i] = np.sum(
                chr_data[int(i * scale) : int(i * scale + scale)]
            )
            return_sample[chr_name] = scaled_chr
    return return_sample


def gender_correct(sample, gender):
    """
    Levels gonosomal reads with the one at the autosomes.
    """

    if gender == "M":
        sample["23"] = sample["23"] * 2
        sample["24"] = sample["24"] * 2

    return sample


def get_z_score(results_c, results):
    """
    Calculates between sample z-score.
    """
    results_nr, results_r, results_w = (
        results["results_nr"],
        results["results_r"],
        results["results_w"],
    )
    zs = []
    for segment in results_c:
        segment_nr = results_nr[segment[0]][segment[1] : segment[2] + 1]
        segment_rr = results_r[segment[0]][segment[1] : segment[2] + 1]
        segment_nr = [
            segment_nr[i] for i in range(len(segment_nr)) if segment_rr[i] != 0
        ]
        for i in range(len(segment_nr)):
            for ii in range(len(segment_nr[i])):
                if not np.isfinite(segment_nr[i][ii]):
                    segment_nr[i][ii] = np.nan
        segment_w = results_w[segment[0]][segment[1] : segment[2] + 1]
        segment_w = [
            segment_w[i] for i in range(len(segment_w)) if segment_rr[i] != 0
        ]
        null_segments = [
            np.ma.average(
                np.ma.masked_array(x, pd.isnull(x)), weights=segment_w
            )
            for x in np.transpose(segment_nr)
        ]
        null_mean = np.ma.mean([x for x in null_segments if np.isfinite(x)])
        null_sd = np.ma.std([x for x in null_segments if np.isfinite(x)])
        z = (segment[3] - null_mean) / null_sd
        z = min(z, 1000)
        z = max(z, -1000)
        if math.isnan(null_mean) or math.isnan(null_sd):
            z = "nan"
        zs.append(z)
    return zs


def get_median_segment_variance(results_c, results_r):
    """
    Returns MSV, measure for sample-wise noise.
    """
    vars = []
    for segment in results_c:
        segment_r = results_r[segment[0]][
            int(segment[1]) : int(segment[2]) + 1
        ]
        segment_r = [x for x in segment_r if x != 0]
        if segment_r:
            var = np.var(segment_r)
            vars.append(var)
    return np.median(vars)
