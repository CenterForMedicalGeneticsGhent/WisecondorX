# WisecondorX

import bisect
import logging
import random
import os
import sys
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
from scipy.signal import argrelextrema
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
import typer

from wisecondorx.overall_tools import scale_sample, gender_correct
from wisecondorx.newref_control import (
    tool_newref_prep,
    tool_newref_main,
    tool_newref_merge,
)

"""
A Gaussian mixture model is fitted against
all one-dimensional reference y-fractions.
Two components are expected: one for males,
and one for females. The local minimum will
serve as the cut-off point.
"""


def wcx_newref(
    infiles: List[Path] = typer.Argument(
        ...,
        help="Path to all reference data files (e.g. path/to/reference/*.npz)",
    ),
    outfile: Path = typer.Argument(
        ...,
        help="Path and filename for the reference output (e.g. path/to/myref.npz)",
    ),
    nipt: bool = typer.Option(
        False, "--nipt", help="Use flag for NIPT", is_flag=True
    ),
    yfrac: float = typer.Option(
        None,
        "--yfrac",
        help="Use to manually set the Y read fraction cutoff, which defines gender",
    ),
    plotyfrac: Path = typer.Option(
        None,
        "--plotyfrac",
        help="Path to yfrac .png plot for optimization; software will stop after plotting",
    ),
    refsize: int = typer.Option(
        300, "--refsize", help="Amount of reference locations per target"
    ),
    binsize: int = int(
        1e5
    ),  # Cannot use scientific notation as default for Typer Option int
    cpus: int = typer.Option(
        1, "--cpus", help="Use multiple cores to find reference bins"
    ),
) -> None:
    """
    Create a new reference using healthy reference samples.
    """
    # Fix the binsize type from typer, because 1e5 is float
    binsize_val = int(binsize)

    logging.info("Creating new reference")

    if yfrac is not None:
        if yfrac < 0 or yfrac > 1:
            logging.critical(
                "Parameter --yfrac should be a positive number lower than or equal to 1"
            )
            sys.exit(1)

    split_path = list(os.path.split(outfile))
    if split_path[-1].endswith(".npz"):
        split_path[-1] = split_path[-1][:-4]

    basepath = os.path.join(split_path[0], split_path[1])
    prepfile = f"{basepath}_prep.npz"
    prepdatafile = f"{basepath}_prep_data.npy"
    partfile = f"{basepath}_part"
    tmpoutfile_A = f"{basepath}.tmp.A.npz"
    tmpoutfile_F = f"{basepath}.tmp.F.npz"
    tmpoutfile_M = f"{basepath}.tmp.M.npz"

    samples: list[dict[str, np.ndarray]] = []
    logging.info("Importing data ...")
    for infile in infiles:
        logging.info(f"Loading: {infile}")
        npzdata = np.load(infile, encoding="latin1", allow_pickle=True)
        sample = npzdata["sample"].item()
        b_size = int(npzdata["binsize"])
        logging.info(f"Binsize: {int(b_size)}")
        samples.append(scale_sample(sample, b_size, binsize_val))

    samples_array = np.array(samples)
    genders, trained_cutoff = train_gender_model(
        samples=samples_array, yfrac=yfrac, plotyfrac=plotyfrac
    )

    if genders.count("F") < 5 and nipt:
        logging.warning(
            "A NIPT reference should have at least 5 female feti samples. Removing --nipt flag."
        )
        nipt = False

    if not nipt:
        for i, sample in enumerate(samples_array):
            samples_array[i] = gender_correct(sample, genders[i])

    total_mask, bins_per_chr = get_mask(samples_array)
    if genders.count("F") > 4:
        mask_F, _ = get_mask(samples_array[np.array(genders) == "F"])
        total_mask = total_mask & mask_F
    if genders.count("M") > 4 and not nipt:
        mask_M, _ = get_mask(samples_array[np.array(genders) == "M"])
        total_mask = total_mask & mask_M

    outfiles_list: List[str] = []

    if len(genders) > 9:
        logging.info("Starting autosomal reference creation ...")
        outfiles_list.append(tmpoutfile_A)
        tool_newref_prep(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            binsize=binsize_val,
            samples=samples_array,
            gender="A",
            mask=total_mask,
            bins_per_chr=bins_per_chr,
        )
        logging.info("This might take a while ...")
        tool_newref_main(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            partfile=partfile,
            tmpoutfile=tmpoutfile_A,
            refsize=refsize,
            cpus=cpus,
        )
    else:
        logging.critical(
            "Provide at least 10 samples to enable the generation of a reference."
        )
        sys.exit(1)

    if genders.count("F") > 4:
        logging.info("Starting female gonosomal reference creation ...")
        outfiles_list.append(tmpoutfile_F)
        tool_newref_prep(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            binsize=binsize_val,
            samples=samples_array[np.array(genders) == "F"],
            gender="F",
            mask=total_mask,
            bins_per_chr=bins_per_chr,
        )
        logging.info("This might take a while ...")
        tool_newref_main(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            partfile=partfile,
            tmpoutfile=tmpoutfile_F,
            refsize=refsize,
            cpus=1,
        )
    else:
        logging.warning(
            "Provide at least 5 female samples to enable normalization of female gonosomes."
        )

    if not nipt:
        if genders.count("M") > 4:
            logging.info("Starting male gonosomal reference creation ...")
            outfiles_list.append(tmpoutfile_M)
            tool_newref_prep(
                prepdatafile=prepdatafile,
                prepfile=prepfile,
                binsize=binsize_val,
                samples=samples_array[np.array(genders) == "M"],
                gender="M",
                mask=total_mask,
                bins_per_chr=bins_per_chr,
            )
            tool_newref_main(
                prepdatafile=prepdatafile,
                prepfile=prepfile,
                partfile=partfile,
                tmpoutfile=tmpoutfile_M,
                refsize=refsize,
                cpus=1,
            )
        else:
            logging.warning(
                "Provide at least 5 male samples to enable normalization of male gonosomes."
            )

    tool_newref_merge(
        outfile=outfile,
        nipt=nipt,
        outfiles=outfiles_list,
        trained_cutoff=trained_cutoff,
    )

    logging.info("Finished creating reference")


def train_gender_model(
    samples: np.ndarray,
    yfrac: Optional[float] = None,
    plotyfrac: Optional[str] = None,
) -> Tuple[List[str], float]:
    genders = np.empty(len(samples), dtype="object")
    y_fractions = []
    for sample in samples:
        y_fractions.append(
            float(np.sum(sample["24"]))
            / float(np.sum([np.sum(sample[x]) for x in sample.keys()]))
        )
    y_fractions = np.array(y_fractions)

    gmm = GaussianMixture(
        n_components=2,
        covariance_type="full",
        reg_covar=1e-99,
        max_iter=10000,
        tol=1e-99,
    )
    gmm.fit(X=y_fractions.reshape(-1, 1))
    gmm_x = np.linspace(0, 0.02, 5000)
    gmm_y = np.exp(gmm.score_samples(gmm_x.reshape(-1, 1)))

    if plotyfrac is not None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(16, 6))
        ax.hist(y_fractions, bins=100, density=True)
        ax.plot(gmm_x, gmm_y, "r-", label="Gaussian mixture fit")
        ax.set_xlim([0, 0.02])
        ax.legend(loc="best")
        plt.savefig(plotyfrac)
        logging.info("Image written to {}, now quitting ...".format(plotyfrac))
        exit()

    if yfrac is not None:
        cut_off = yfrac
    else:
        sort_idd = np.argsort(gmm_x)
        sorted_gmm_y = gmm_y[sort_idd]

        local_min_i = argrelextrema(sorted_gmm_y, np.less)

        cut_off = gmm_x[local_min_i][0]
        logging.info(
            "Determined --yfrac cutoff: {}".format(str(round(cut_off, 4)))
        )

    genders[y_fractions > cut_off] = "M"
    genders[y_fractions < cut_off] = "F"

    return genders.tolist(), cut_off


"""
Finds mask (locations of bins without data) in the
subset 'samples'.
"""


def get_mask(samples: np.ndarray) -> Tuple[np.ndarray, List[int]]:
    by_chr = []
    bins_per_chr = []
    sample_count = len(samples)

    for chr in range(1, 25):
        max_len = max([sample[str(chr)].shape[0] for sample in samples])
        this_chr = np.zeros((max_len, sample_count), dtype=float)
        bins_per_chr.append(max_len)
        i = 0
        for sample in samples:
            this_chr[:, i] = sample[str(chr)]
            i += 1
        by_chr.append(this_chr)
    all_data = np.concatenate(by_chr, axis=0)

    sum_per_sample = np.sum(all_data, 0)
    all_data = all_data / sum_per_sample

    sum_per_bin = np.sum(all_data, 1)
    mask = sum_per_bin > 0

    return mask, bins_per_chr


"""
Normalizes samples for read depth and applies mask.
"""


def normalize_and_mask(
    samples: np.ndarray, chrs: List[int], mask: np.ndarray
) -> np.ndarray:
    by_chr = []
    sample_count = len(samples)

    for chr in chrs:
        max_len = max([sample[str(chr)].shape[0] for sample in samples])
        this_chr = np.zeros((max_len, sample_count), dtype=float)
        i = 0
        for sample in samples:
            this_chr[:, i] = sample[str(chr)]
            i += 1
        by_chr.append(this_chr)
    all_data = np.concatenate(by_chr, axis=0)

    sum_per_sample = np.sum(all_data, 0)
    all_data = all_data / sum_per_sample

    masked_data = all_data[mask, :]

    return masked_data


"""
Executes PCA. Rotations are saved which enable
between sample normalization in the test phase.
"""


def train_pca(
    ref_data: np.ndarray, pcacomp: int = 5
) -> Tuple[np.ndarray, PCA]:
    t_data = ref_data.T
    pca = PCA(n_components=pcacomp)
    pca.fit(t_data)
    PCA(copy=True, whiten=False)
    transformed = pca.transform(t_data)
    inversed = pca.inverse_transform(transformed)
    corrected = t_data / inversed

    return corrected.T, pca


"""
Calculates within-sample reference.
"""


def get_reference(
    pca_corrected_data: np.ndarray,
    masked_bins_per_chr: List[int],
    masked_bins_per_chr_cum: List[int],
    ref_size: int,
    part: int,
    split_parts: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    big_indexes = []
    big_distances = []

    bincount = masked_bins_per_chr_cum[-1]

    start_num, end_num = _get_part(part - 1, split_parts, bincount)
    logging.info(
        "Working on thread {} of {}, meaning bins {} up to {}".format(
            part, split_parts, start_num, end_num
        )
    )
    regions = _split_by_chr(start_num, end_num, masked_bins_per_chr_cum)

    for region in regions:
        chr = region[0]
        start = region[1]
        end = region[2]

        if start_num > start:
            start = start_num
        if end_num < end:
            end = end_num

        if len(masked_bins_per_chr_cum) > 22 and chr != 22 and chr != 23:
            part_indexes = np.zeros((end - start, ref_size), dtype=np.int32)
            part_distances = np.ones((end - start, ref_size))
            big_indexes.extend(part_indexes)
            big_distances.extend(part_distances)
            continue
        chr_data = np.concatenate(
            (
                pca_corrected_data[
                    : masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr],
                    :,
                ],
                pca_corrected_data[masked_bins_per_chr_cum[chr] :, :],
            )
        )

        part_indexes, part_distances = get_ref_for_bins(
            ref_size, start, end, pca_corrected_data, chr_data
        )

        big_indexes.extend(part_indexes)
        big_distances.extend(part_distances)

    index_array = np.array(big_indexes)
    distance_array = np.array(big_distances)
    null_ratio_array = np.zeros(
        (len(distance_array), min(len(pca_corrected_data[0]), 100))
    )
    samples = np.transpose(pca_corrected_data)
    for null_i, case_i in enumerate(
        random.sample(
            range(len(pca_corrected_data[0])),
            min(len(pca_corrected_data[0]), 100),
        )
    ):
        sample = samples[case_i]
        for bin_i in list(range(len(sample)))[start_num:end_num]:
            ref = sample[index_array[bin_i - start_num]]
            r = np.log2(sample[bin_i] / np.median(ref))
            null_ratio_array[bin_i - start_num][null_i] = r
    return index_array, distance_array, null_ratio_array


def _split_by_chr(
    start: int, end: int, chr_bin_sums: List[int]
) -> List[List[int]]:
    areas = []
    tmp = [0, start, 0]
    for i, val in enumerate(chr_bin_sums):
        tmp[0] = i
        if val >= end:
            break
        if start < val < end:
            tmp[2] = val
            areas.append(tmp)
            tmp = [i, val, 0]
        tmp[1] = val
    tmp[2] = end
    areas.append(tmp)
    return areas


def _get_part(partnum: int, outof: int, bincount: int) -> Tuple[int, int]:
    start_bin = int(bincount / float(outof) * partnum)
    end_bin = int(bincount / float(outof) * (partnum + 1))
    return start_bin, end_bin


"""
Calculates within-sample reference for a particular chromosome.
"""


def get_ref_for_bins(
    ref_size: int,
    start: int,
    end: int,
    pca_corrected_data: np.ndarray,
    chr_data: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    find_pos = bisect.bisect
    ref_indexes = np.zeros((end - start, ref_size), dtype=np.int32)
    ref_distances = np.ones((end - start, ref_size))
    for this_bin in range(start, end):
        this_mask = np.sum(
            np.power(chr_data - pca_corrected_data[this_bin, :], 2), 1
        )
        this_indexes = [-1 for i in range(ref_size)]
        this_distances = [1e10 for i in range(ref_size)]
        remove_index = this_indexes.pop
        remove_dist = this_distances.pop
        insert_index = this_indexes.insert
        insert_dist = this_distances.insert
        cur_max = 1e10
        for i, binVal in enumerate(this_mask):
            if binVal < cur_max:
                pos = find_pos(this_distances, binVal)
                remove_index(-1)
                remove_dist(-1)
                insert_index(pos, i)
                insert_dist(pos, binVal)
                cur_max = this_distances[-1]
        ref_indexes[this_bin - start, :] = this_indexes
        ref_distances[this_bin - start, :] = this_distances
    return ref_indexes, ref_distances
