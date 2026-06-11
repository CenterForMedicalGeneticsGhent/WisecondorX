# WisecondorX

import copy
import logging
import os
import sys
import time

import numpy as np
from concurrent import futures
import bisect
import random
from pathlib import Path
import typer
from scipy.signal import argrelextrema
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from wisecondorx.ref_qc import qc_reference
from wisecondorx.overall_tools import gender_correct, scale_sample
import argparse
from typing import Annotated


"""
Outputs preparation files of read depth normalized
data and contains PCA information to execute between-
sample normalization during testing. Function is
executed three times. Once for autosomes, once for XX
gonosomes (if enough females are included) and once
for XY gonosomes (if enough males are included).
"""


def _robust_cutoff(distances, mad_multiplier, floor):
    """
    Computes a robust distance cutoff to remove unmatchable bins.

    :param mad_multiplier: Multiplier applied to the median absolute deviation (MAD).
    :param floor: Absolute minimum value for the cutoff.
    """
    med = np.median(distances)
    mad = np.median(np.abs(distances - med))
    return max(med + mad_multiplier * mad, floor)


def _masked_chr_ids(mask, bins_per_chr):
    """
    Maps each `True` index in the 1D mask back to its 1-based chromosome ID.
    Used to split PCA distance filtering thresholds by chromosome class.
    """
    masked_indices = np.where(mask)[0]
    chr_ends = np.cumsum(bins_per_chr)
    chr_ids = np.searchsorted(chr_ends, masked_indices, side="right") + 1
    return masked_indices, chr_ids


def tool_newref_prep(args, samples, gender, mask, bins_per_chr):
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

    masked_indices, masked_chr_ids = _masked_chr_ids(mask, bins_per_chr)
    bad_bins_mask = np.zeros(len(masked_indices), dtype=bool)

    # PCA Distance filtering: calculate outlier thresholds independently per chromosome class.
    # This prevents the noisier chrX/chrY from being heavily pruned when judged against an
    # autosome-dominated median profile.
    #
    # We use higher MAD multipliers (20.0-50.0) and floors (10.0-15.0)
    # than typical outlier detection to ensure the reference retains
    # enough natural variance to prevent massive NaN gaps when `predict` runs on shallower,
    # noisier routine samples.
    if gender == "A":
        class_settings = [("autosomes", masked_chr_ids <= 22, 20.0, 10.0)]
    elif gender == "F":
        class_settings = [("chrX", masked_chr_ids == 23, 20.0, 10.0)]
    else:
        class_settings = [
            ("chrX", masked_chr_ids == 23, 20.0, 10.0),
            ("chrY", masked_chr_ids == 24, 50.0, 15.0),
        ]
    for class_name, class_mask, mad_mult, floor in class_settings:
        if np.sum(class_mask) < 10:
            continue
        class_data = pca_corrected_data[class_mask]
        med_prof = np.median(class_data, axis=0)
        dist_to_med = np.sum((class_data - med_prof) ** 2, axis=1)
        cutoff = _robust_cutoff(dist_to_med, mad_mult, floor)
        class_bad_bins = dist_to_med > cutoff
        bad_bins_mask[class_mask] = class_bad_bins
        if np.any(class_bad_bins):
            logging.info(
                "Removing {} anomalous {} bins based on PCA distance (cutoff={:.4f})".format(
                    int(np.sum(class_bad_bins)), class_name, cutoff
                )
            )

    if np.any(bad_bins_mask):
        global_bad_indices = masked_indices[bad_bins_mask]
        mask[global_bad_indices] = False

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

    np.save(args.prepdatafile, pca_corrected_data)

    np.savez_compressed(
        args.prepfile,
        binsize=args.binsize,
        gender=gender,
        mask=mask,
        bins_per_chr=bins_per_chr,
        masked_bins_per_chr=masked_bins_per_chr,
        masked_bins_per_chr_cum=masked_bins_per_chr_cum,
        pca_components=pca.components_,
        pca_mean=pca.mean_,
    )


def tool_newref_main(args, cpus):
    """
    Prepares subfiles if multi-threading is requested.
    Main file is split in 'cpus' subfiles, each subfile
    is processed by a separate thread.
    """
    pca_corrected_data = np.load(args.prepdatafile, mmap_mode="r")
    if cpus != 1:
        with futures.ThreadPoolExecutor(max_workers=args.cpus) as executor:
            for part in range(1, cpus + 1):
                this_args = copy.copy(args)
                this_args.part = [part, cpus]
                executor.submit(
                    _tool_newref_part, this_args, pca_corrected_data
                )
            executor.shutdown(wait=True)
    else:
        for part in range(1, cpus + 1):
            args.part = [part, cpus]
            _tool_newref_part(args, pca_corrected_data)

    tool_newref_post(args, cpus)

    os.remove(args.prepfile)
    os.remove(args.prepdatafile)
    for part in range(1, cpus + 1):
        os.remove("{}_{}.npz".format(args.partfile, str(part)))


def _tool_newref_part(args, pca_corrected_data):
    """
    Function executed once for each thread. Controls
    within-sample reference creation.
    """
    if args.part[0] > args.part[1]:
        logging.critical(
            "Part should be smaller or equal to total parts:{} > {} is wrong".format(
                args.part[0], args.part[1]
            )
        )
        sys.exit()
    if args.part[0] < 0:
        logging.critical(
            "Part should be at least zero: {} < 0 is wrong".format(
                args.part[0]
            )
        )
        sys.exit()

    npzdata = np.load(args.prepfile, encoding="latin1", allow_pickle=True)
    masked_bins_per_chr = npzdata["masked_bins_per_chr"]
    masked_bins_per_chr_cum = npzdata["masked_bins_per_chr_cum"]

    indexes, distances, null_ratios = get_reference(
        pca_corrected_data,
        masked_bins_per_chr,
        masked_bins_per_chr_cum,
        ref_size=args.refsize,
        part=args.part[0],
        split_parts=args.part[1],
    )

    np.savez_compressed(
        "{}_{}.npz".format(args.partfile, str(args.part[0])),
        indexes=indexes,
        distances=distances,
        null_ratios=null_ratios,
    )


def tool_newref_post(args, cpus):
    """
    Merges separate subfiles (one for each thread) to a
    new temporary output file.
    """
    npzdata_prep = np.load(args.prepfile, encoding="latin1", allow_pickle=True)

    big_indexes = []
    big_distances = []
    big_null_ratios = []
    for part in range(1, cpus + 1):
        infile = "{}_{}.npz".format(args.partfile, str(part))
        npzdata_part = np.load(infile, encoding="latin1")
        big_indexes.extend(npzdata_part["indexes"])
        big_distances.extend(npzdata_part["distances"])
        big_null_ratios.extend(npzdata_part["null_ratios"])

    indexes = np.array(big_indexes)
    distances = np.array(big_distances)
    null_ratios = np.array(big_null_ratios)

    np.savez_compressed(
        args.tmpoutfile,
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


def force_remove(file_id):
    """
    Tries to remove text file, when it is busy, until becomes successful.
    This function, prevents OSError: [Errno 26] Text file busy...
    """
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


def tool_newref_merge(args, outfiles, trained_cutoff):
    """
    Merges separate subfiles (A, F, M) to one final
    reference file.
    """
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
    final_ref["is_nipt"] = args.nipt
    final_ref["trained_cutoff"] = trained_cutoff
    np.savez_compressed(args.outfile, **final_ref)


#### newref_tools


def train_gender_model(args, samples):
    """
    A Gaussian mixture model is fitted against
    all one-dimensional reference y-fractions.
    Two components are expected: one for males,
    and one for females. The local minimum will
    serve as the cut-off point.
    """

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

    if args.plotyfrac is not None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(16, 6))
        ax.hist(y_fractions, bins=100, density=True)
        ax.plot(gmm_x, gmm_y, "r-", label="Gaussian mixture fit")
        ax.set_xlim([0, 0.02])
        ax.legend(loc="best")
        plt.savefig(args.plotyfrac)
        logging.info(
            "Image written to {}, now quitting ...".format(args.plotyfrac)
        )
        sys.exit()

    if args.yfrac is not None:
        cut_off = args.yfrac
    else:
        sort_idd = np.argsort(gmm_x)
        sorted_gmm_y = gmm_y[sort_idd]

        local_min_i = argrelextrema(sorted_gmm_y, np.less)

        if len(local_min_i[0]) > 0:
            cut_off = gmm_x[local_min_i][0]
        else:
            cut_off = float(np.mean(gmm.means_))
            logging.warning(
                "No local minimum found, using mean of GMM components as fallback"
            )

        logging.info(
            "Determined --yfrac cutoff: {}".format(str(round(cut_off, 4)))
        )

    genders[y_fractions > cut_off] = "M"
    genders[y_fractions < cut_off] = "F"

    return genders.tolist(), cut_off


def get_mask(samples):
    """
    Finds mask (locations of bins without data) in the
    subset 'samples'.
    """
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
    # Mask out bins with 0 reads OR reads less than 5% of the median coverage
    # to reduce noise from small binsizes
    median_cov = np.median(sum_per_bin[sum_per_bin > 0])
    mask = sum_per_bin > (0.05 * median_cov)

    return mask, bins_per_chr


def normalize_and_mask(samples, chrs, mask):
    """
    Normalizes samples for read depth and applies mask.
    """
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


def train_pca(ref_data, pcacomp=5):
    """
    Executes PCA. Rotations are saved which enable
    between sample normalization in the test phase.
    """
    t_data = ref_data.T
    pca = PCA(n_components=pcacomp)
    pca.fit(t_data)
    transformed = pca.transform(t_data)
    inversed = pca.inverse_transform(transformed)
    corrected = t_data / inversed

    return corrected.T, pca


def get_reference(
    pca_corrected_data,
    masked_bins_per_chr,
    masked_bins_per_chr_cum,
    ref_size,
    part,
    split_parts,
):
    """
    Calculates within-sample reference.
    """
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


def _split_by_chr(start, end, chr_bin_sums):
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


def _get_part(partnum, outof, bincount):
    start_bin = int(bincount / float(outof) * partnum)
    end_bin = int(bincount / float(outof) * (partnum + 1))
    return start_bin, end_bin


def get_ref_for_bins(ref_size, start, end, pca_corrected_data, chr_data):
    """
    Calculates within-sample reference for a particular chromosome.
    """
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


### main function
def _build_branch_masks(samples, genders, is_nipt):
    total_mask, bins_per_chr = get_mask(samples)

    bins_per_chr = np.array(bins_per_chr, dtype=int)
    chr22_end = int(np.sum(bins_per_chr[:22]))
    chr23_end = int(np.sum(bins_per_chr[:23]))
    chr24_end = int(np.sum(bins_per_chr[:24]))

    mask_a = np.zeros_like(total_mask, dtype=bool)
    mask_a[:chr22_end] = total_mask[:chr22_end]

    mask_f_sex = None
    if genders.count("F") > 4:
        female_mask, _ = get_mask(samples[np.array(genders) == "F"])
        mask_f_sex = np.zeros_like(total_mask, dtype=bool)
        if len(female_mask) >= chr23_end:
            mask_f_sex[chr22_end:chr23_end] = female_mask[chr22_end:chr23_end]

    mask_m_sex = None
    if genders.count("M") > 4 and not is_nipt:
        male_mask, _ = get_mask(samples[np.array(genders) == "M"])
        mask_m_sex = np.zeros_like(total_mask, dtype=bool)
        if len(male_mask) >= chr23_end:
            mask_m_sex[chr22_end:chr23_end] = male_mask[chr22_end:chr23_end]
        if len(male_mask) >= chr24_end:
            mask_m_sex[chr23_end:chr24_end] = male_mask[chr23_end:chr24_end]

    return bins_per_chr.tolist(), mask_a, mask_f_sex, mask_m_sex


def wcx_newref(
    infiles: list[Path] = typer.Argument(
        ...,
        help="Path to all reference data files (e.g. path/to/reference/*.npz)",
    ),
    outfile: Path = typer.Argument(
        ...,
        help="Path and filename for the reference output (e.g. path/to/myref.npz)",
    ),
    nipt: bool = typer.Option(False, "--nipt", help="Use flag for NIPT"),
    yfrac: Annotated[
        float,
        typer.Option(
            "--yfrac",
            min=0.0,
            max=1.0,
            help="Use to manually set the Y read fraction cutoff, which defines sex",
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
    binsize: int = typer.Option(
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

    args = argparse.Namespace()
    args.infiles = infiles
    args.outfile = outfile
    args.nipt = nipt
    args.yfrac = yfrac
    args.plotyfrac = plotyfrac
    args.refsize = refsize
    args.binsize = binsize
    args.cpus = cpus

    if args.yfrac is not None:
        if args.yfrac < 0 or args.yfrac > 1:
            logging.critical(
                "Parameter --yfrac should be a positive number lower than or equal to 1"
            )
            sys.exit()

    split_path = list(os.path.split(args.outfile))
    if split_path[-1][-4:] == ".npz":
        split_path[-1] = split_path[-1][:-4]
    base_path = os.path.join(split_path[0], split_path[1])

    args.basepath = base_path
    args.prepfile = "{}_prep.npz".format(base_path)
    args.prepdatafile = "{}_prep_data.npy".format(base_path)
    args.partfile = "{}_part".format(base_path)

    samples = []
    logging.info("Importing data ...")
    for infile in args.infiles:
        logging.info("Loading: {}".format(infile))
        npzdata = np.load(infile, encoding="latin1", allow_pickle=True)
        sample = npzdata["sample"].item()
        binsize = int(npzdata["binsize"])
        logging.info("Binsize: {}".format(int(binsize)))
        samples.append(scale_sample(sample, binsize, args.binsize))

    samples = np.array(samples)
    genders, trained_cutoff = train_gender_model(args, samples)

    if genders.count("F") < 5 and args.nipt:
        logging.warning(
            "A NIPT reference should have at least 5 female feti samples. Removing --nipt flag."
        )
        args.nipt = False
    if not args.nipt:
        for i, sample in enumerate(samples):
            samples[i] = gender_correct(sample, genders[i])

    bins_per_chr, mask_A, mask_F_sex, mask_M_sex = _build_branch_masks(
        samples, genders, args.nipt
    )
    chr22_end = int(np.sum(bins_per_chr[:22]))
    chr23_end = int(np.sum(bins_per_chr[:23]))
    chr24_end = int(np.sum(bins_per_chr[:24]))

    outfiles = []
    if len(genders) > 9:
        logging.info("Starting autosomal reference creation ...")
        args.tmpoutfile = "{}.tmp.A.npz".format(args.basepath)
        outfiles.append(args.tmpoutfile)
        tool_newref_prep(args, samples, "A", mask_A, bins_per_chr)
        logging.info("This might take a while ...")
        tool_newref_main(args, args.cpus)
    else:
        logging.critical(
            "Provide at least 10 samples to enable the generation of a reference."
        )
        sys.exit()

    if genders.count("F") > 4:
        logging.info("Starting female gonosomal reference creation ...")
        args.tmpoutfile = "{}.tmp.F.npz".format(args.basepath)
        outfiles.append(args.tmpoutfile)
        mask_F = np.copy(mask_A)
        if mask_F_sex is not None:
            mask_F[chr22_end:chr23_end] = mask_F_sex[chr22_end:chr23_end]
        tool_newref_prep(
            args, samples[np.array(genders) == "F"], "F", mask_F, bins_per_chr
        )
        logging.info("This might take a while ...")
        tool_newref_main(args, 1)
    else:
        logging.warning(
            "Provide at least 5 female samples to enable normalization of female gonosomes."
        )

    if not args.nipt:
        if genders.count("M") > 4:
            logging.info("Starting male gonosomal reference creation ...")
            args.tmpoutfile = "{}.tmp.M.npz".format(args.basepath)
            outfiles.append(args.tmpoutfile)
            mask_M = np.copy(mask_A)
            if mask_M_sex is not None:
                mask_M[chr22_end:chr24_end] = mask_M_sex[chr22_end:chr24_end]
            tool_newref_prep(
                args,
                samples[np.array(genders) == "M"],
                "M",
                mask_M,
                bins_per_chr,
            )
            tool_newref_main(args, 1)
        else:
            logging.warning(
                "Provide at least 5 male samples to enable normalization of male gonosomes."
            )

    tool_newref_merge(args, outfiles, trained_cutoff)

    logging.info("Running QC on the newly created reference...")
    qc_reference(args.outfile)

    logging.info("Finished creating reference")
