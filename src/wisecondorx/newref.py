# WisecondorX

import bisect
import logging
import random
import os
import sys
import time
from pathlib import Path
from typing import Optional, Annotated

import numpy as np
from concurrent import futures
from scipy.signal import argrelextrema
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
import typer

from wisecondorx.utils import scale_sample, gender_correct


def tool_newref_prep(
    prepdatafile: str,
    prepfile: str,
    binsize: int,
    samples: np.ndarray,
    gender: str,
    mask: np.ndarray,
    bins_per_chr: list[int],
) -> None:
    from wisecondorx.newref_tools import (
        normalize_and_mask,
        train_pca,
    )

    if gender == "A":
        last_chr = 22
    elif gender == "F":
        last_chr = 23
    else:
        last_chr = 24

    bins_per_chr = bins_per_chr[:last_chr]
    mask = mask[: np.sum(bins_per_chr)]

    masked_data = normalize_and_mask(samples, range(1, last_chr + 1), mask)
    pca_corrected_data, pca = train_pca(masked_data)

    masked_bins_per_chr = [
        sum(mask[sum(bins_per_chr[:i]) : sum(bins_per_chr[:i]) + x])
        for i, x in enumerate(bins_per_chr)
    ]
    masked_bins_per_chr_cum = [
        sum(masked_bins_per_chr[: x + 1])
        for x in range(len(masked_bins_per_chr))
    ]

    np.save(prepdatafile, pca_corrected_data)

    np.savez_compressed(
        prepfile,
        binsize=binsize,
        gender=gender,
        mask=mask,
        bins_per_chr=bins_per_chr,
        masked_bins_per_chr=masked_bins_per_chr,
        masked_bins_per_chr_cum=masked_bins_per_chr_cum,
        pca_components=pca.components_,
        pca_mean=pca.mean_,
    )


"""
Prepares subfiles if multi-threading is requested.
Main file is split in 'cpus' subfiles, each subfile
is processed by a separate thread.
"""


def tool_newref_main(
    prepdatafile: str,
    prepfile: str,
    partfile: str,
    tmpoutfile: str,
    refsize: int,
    cpus: int,
    part: Optional[list[int]] = None,
) -> None:
    pca_corrected_data = np.load(prepdatafile, mmap_mode="r")
    if cpus != 1:
        with futures.ThreadPoolExecutor(max_workers=cpus) as executor:
            for p in range(1, cpus + 1):
                executor.submit(
                    _tool_newref_part,
                    prepfile,
                    partfile,
                    refsize,
                    [p, cpus],
                    pca_corrected_data,
                )
            executor.shutdown(wait=True)
    else:
        for p in range(1, cpus + 1):
            _tool_newref_part(
                prepfile, partfile, refsize, [p, cpus], pca_corrected_data
            )

    tool_newref_post(prepfile, partfile, tmpoutfile, cpus)

    os.remove(prepfile)
    os.remove(prepdatafile)
    for p in range(1, cpus + 1):
        os.remove("{}_{}.npz".format(partfile, str(p)))


"""
Function executed once for each thread. Controls
within-sample reference creation.
"""


def _tool_newref_part(
    prepfile: str,
    partfile: str,
    refsize: int,
    part: list[int],
    pca_corrected_data: np.ndarray,
) -> None:
    from wisecondorx.newref_tools import get_reference

    if part[0] > part[1]:
        logging.critical(
            "Part should be smaller or equal to total parts:{} > {} is wrong".format(
                part[0], part[1]
            )
        )
        sys.exit()
    if part[0] < 0:
        logging.critical(
            "Part should be at least zero: {} < 0 is wrong".format(part[0])
        )
        sys.exit()

    npzdata = np.load(prepfile, encoding="latin1", allow_pickle=True)
    masked_bins_per_chr = npzdata["masked_bins_per_chr"]
    masked_bins_per_chr_cum = npzdata["masked_bins_per_chr_cum"]

    indexes, distances, null_ratios = get_reference(
        pca_corrected_data,
        masked_bins_per_chr,
        masked_bins_per_chr_cum,
        ref_size=refsize,
        part=part[0],
        split_parts=part[1],
    )

    np.savez_compressed(
        "{}_{}.npz".format(partfile, str(part[0])),
        indexes=indexes,
        distances=distances,
        null_ratios=null_ratios,
    )


"""
Merges separate subfiles (one for each thread) to a
new temporary output file.
"""


def tool_newref_post(
    prepfile: str, partfile: str, tmpoutfile: str, cpus: int
) -> None:
    npzdata_prep = np.load(prepfile, encoding="latin1", allow_pickle=True)

    big_indexes = []
    big_distances = []
    big_null_ratios = []
    for p in range(1, cpus + 1):
        infile = "{}_{}.npz".format(partfile, str(p))
        npzdata_part = np.load(infile, encoding="latin1")
        big_indexes.extend(npzdata_part["indexes"])
        big_distances.extend(npzdata_part["distances"])
        big_null_ratios.extend(npzdata_part["null_ratios"])

    indexes = np.array(big_indexes)
    distances = np.array(big_distances)
    null_ratios = np.array(big_null_ratios)

    np.savez_compressed(
        tmpoutfile,
        binsize=npzdata_prep["binsize"].item(),
        gender=npzdata_prep["gender"].item(),
        mask=npzdata_prep["mask"],
        bins_per_chr=npzdata_prep["bins_per_chr"],
        masked_bins_per_chr=npzdata_prep["masked_bins_per_chr"],
        masked_bins_per_chr_cum=npzdata_prep["masked_bins_per_chr_cum"],
        pca_components=npzdata_prep["pca_components"],
        pca_mean=npzdata_prep["pca_mean"],
        indexes=indexes,
        distances=distances,
        null_ratios=null_ratios,
    )


"""
Tries to remove text file, when it is busy, until becomes successful.
This function, prevents OSError: [Errno 26] Text file busy...
"""


def force_remove(file_id: str) -> None:
    attemp = 1
    while True:
        try:
            os.remove(file_id)
            break
        except Exception:
            print(
                "Attemp #{}: Cannot remove {}, because it is busy, trying again...".format(
                    attemp, file_id
                )
            )
            attemp = attemp + 1
            time.sleep(5)


"""
Merges separate subfiles (A, F, M) to one final
reference file.
"""


def tool_newref_merge(
    outfile: str, nipt: bool, outfiles: list[str], trained_cutoff: float
) -> None:
    final_ref = {"has_female": False, "has_male": False}
    for file_id in outfiles:
        npz_file = np.load(file_id, encoding="latin1", allow_pickle=True)
        gender = str(npz_file["gender"])
        for component in [x for x in npz_file.keys() if x != "gender"]:
            if gender == "F":
                final_ref["has_female"] = True
                final_ref["{}.F".format(str(component))] = npz_file[component]
            elif gender == "M":
                final_ref["has_male"] = True
                final_ref["{}.M".format(str(component))] = npz_file[component]
            else:
                final_ref[str(component)] = npz_file[component]
        force_remove(file_id)
    final_ref["is_nipt"] = nipt
    final_ref["trained_cutoff"] = trained_cutoff
    np.savez_compressed(outfile, **final_ref)


"""
A Gaussian mixture model is fitted against
all one-dimensional reference y-fractions.
Two components are expected: one for males,
and one for females. The local minimum will
serve as the cut-off point.
"""


def wcx_newref(
    infiles: list[Path] = typer.Argument(
        ...,
        help="Path to all reference data files (e.g. path/to/reference/*.npz)",
    ),
    prefix: Path = typer.Argument(
        ...,
        help="Prefix for the reference output (e.g. path/to/myref)",
    ),
    nipt: bool = typer.Option(False, "--nipt", help="Use flag for NIPT"),
    yfrac: Annotated[
        float,
        typer.Option(
            min=0.0,
            max=1.0,
            help="Use to manually set the Y read fraction cutoff, which defines gender",
        ),
    ] = None,
    plotyfrac: Path = typer.Option(
        None,
        "--plotyfrac",
        help="Path to yfrac .png plot for optimization; software will stop after plotting",
    ),
    refsize: int = typer.Option(
        300, "--refsize", help="Amount of reference locations per target"
    ),
    target_binsize: int = typer.Option(
        5000, "--binsize", help="Size of target bins in base pairs"
    ),
    cpus: int = typer.Option(
        1, "--cpus", help="Use multiple cores to find reference bins"
    ),
) -> None:
    """
    Create a new reference using healthy reference samples.
    """

    logging.info("Creating new reference")

    outfile = Path(prefix, ".npz")
    prepfile = Path(prefix, "_prep.npz")
    prepdatafile = Path(prefix, "_prep_data.npy")
    partfile = Path(prefix, "_part")
    tmpoutfile_A = Path(prefix, ".tmp.A.npz")
    tmpoutfile_F = Path(prefix, ".tmp.F.npz")
    tmpoutfile_M = Path(prefix, ".tmp.M.npz")

    samples: list[dict[str, np.ndarray]] = []
    logging.info("Importing data ...")
    for infile in infiles:
        logging.info(f"Loading: {infile}")
        npzdata = np.load(infile, encoding="latin1", allow_pickle=True)
        sample = npzdata["sample"].item()
        source_binsize = int(npzdata["binsize"])
        samples.append(scale_sample(sample, source_binsize, target_binsize))

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

    outfiles_list: list[str] = []

    if len(genders) > 9:
        logging.info("Starting autosomal reference creation ...")
        outfiles_list.append(tmpoutfile_A)
        tool_newref_prep(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            binsize=target_binsize,
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
            binsize=target_binsize,
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
                binsize=target_binsize,
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
    yfrac: float = None,
    plotyfrac: str = None,
) -> tuple[list[str], float]:
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


def get_mask(samples: np.ndarray) -> tuple[np.ndarray, list[int]]:
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
    samples: np.ndarray, chrs: list[int], mask: np.ndarray
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
) -> tuple[np.ndarray, PCA]:
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
    masked_bins_per_chr: list[int],
    masked_bins_per_chr_cum: list[int],
    ref_size: int,
    part: int,
    split_parts: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    start: int, end: int, chr_bin_sums: list[int]
) -> list[list[int]]:
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


def _get_part(partnum: int, outof: int, bincount: int) -> tuple[int, int]:
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
) -> tuple[np.ndarray, np.ndarray]:
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
