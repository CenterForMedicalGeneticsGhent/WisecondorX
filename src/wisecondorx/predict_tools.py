# WisecondorX

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
from wisecondorx.overall_tools import Sex


def predict_gender(sample, trained_cutoff):
    """
    Returns gender based on Gaussian mixture
    model trained during newref phase.
    """
    Y_fraction = float(np.sum(sample["24"])) / float(
        np.sum([np.sum(sample[x]) for x in sample.keys()])
    )
    if Y_fraction > trained_cutoff:
        return Sex.MALE
    else:
        return Sex.FEMALE


"""
Normalize sample for read depth and apply mask.
"""


def coverage_normalize_and_mask(sample, ref_file, ap):
    by_chr = []

    chrs = range(1, len(ref_file["bins_per_chr{}".format(ap)]) + 1)

    for chr in chrs:
        this_chr = np.zeros(
            ref_file["bins_per_chr{}".format(ap)][chr - 1], dtype=float
        )
        min_len = min(
            ref_file["bins_per_chr{}".format(ap)][chr - 1],
            len(sample[str(chr)]),
        )
        this_chr[:min_len] = sample[str(chr)][:min_len]
        by_chr.append(this_chr)
    all_data = np.concatenate(by_chr, axis=0)
    all_data = all_data / np.sum(all_data)
    masked_data = all_data[ref_file["mask{}".format(ap)]]

    return masked_data


"""
Project test sample to PCA space.
"""


def project_pc(sample_data, ref_file, ap):
    components = ref_file["pca_components{}".format(ap)]
    mean = ref_file["pca_mean{}".format(ap)]

    centered_data = sample_data - mean
    transform = np.dot(centered_data, components.T)

    reconstructed = np.dot(transform, components) + mean
    return sample_data / reconstructed


"""
Defines cutoff that will add bins to a blacklist
depending on the within reference distances.
"""


def get_optimal_cutoff(ref_file, repeats):
    distances = ref_file["distances"]
    cutoff = float("inf")
    for i in range(0, repeats):
        mask = distances < cutoff
        average = np.average(distances[mask])
        stddev = np.std(distances[mask])
        cutoff = average + 3 * stddev
    return cutoff


"""
Within sample normalization. Cycles through a number
of repeats where z-scores define whether a bin is seen
as 'normal' in a sample (in most cases this means
'non-aberrant') -- if it is, it can be used as a
reference for the other bins.
"""


def normalize_repeat(test_data, ref_file, optimal_cutoff, ct, cp, ap):
    results_z = None
    results_r = None
    ref_sizes = None
    test_copy = np.copy(test_data)
    for i in range(3):
        results_z, results_r, ref_sizes = _normalize_once(
            test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap
        )

        test_copy[ct:][np.abs(results_z) >= norm.ppf(0.99)] = -1

    if results_r is None or results_z is None or ref_sizes is None:
        raise RuntimeError("Normalization failed to produce output arrays")

    m_lr = np.nanmedian(np.log2(results_r))
    m_z = np.nanmedian(results_z)

    return results_z, results_r, ref_sizes, m_lr, m_z


def _normalize_once(
    test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap
):
    masked_bins_per_chr = ref_file["masked_bins_per_chr{}".format(ap)]
    masked_bins_per_chr_cum = ref_file["masked_bins_per_chr_cum{}".format(ap)]
    results_z = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    results_r = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    ref_sizes = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    indexes = ref_file["indexes{}".format(ap)]
    distances = ref_file["distances{}".format(ap)]

    i = ct
    i2 = 0
    for chr in list(range(len(masked_bins_per_chr)))[cp:]:
        start = masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]
        end = masked_bins_per_chr_cum[chr]
        chr_data = np.concatenate(
            (
                test_copy[
                    : masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]
                ],
                test_copy[masked_bins_per_chr_cum[chr] :],
            )
        )

        for index in indexes[start:end]:
            ref_data = chr_data[index[distances[i] < optimal_cutoff]]
            ref_data = ref_data[ref_data >= 0]
            ref_stdev = np.std(ref_data)
            results_z[i2] = (test_data[i] - np.mean(ref_data)) / ref_stdev
            results_r[i2] = test_data[i] / np.median(ref_data)
            ref_sizes[i2] = ref_data.shape[0]
            i += 1
            i2 += 1

    return results_z, results_r, ref_sizes


"""
The means of sets of within-sample reference
distances can serve as inverse weights for
CBS, Z-scoring and plotting.
"""


def get_weights(ref_file, ap):
    inverse_weights = [
        np.mean(np.sqrt(x)) for x in ref_file["distances{}".format(ap)]
    ]
    weights = np.array([1 / x for x in inverse_weights])
    return weights


"""
Unmasks results array.
"""


def inflate_results(results, rem_input):
    temp = [0 for x in rem_input["mask"]]
    j = 0
    for i, val in enumerate(rem_input["mask"]):
        if val:
            temp[i] = results[j]
            j += 1
    return temp


"""
Log2-transforms results_r. If resulting elements are infinite,
all corresponding possible positions (at results_r, results_z
and results_w are set to 0 (blacklist)).
"""


def log_trans(results, log_r_median):
    for chr in range(len(results["results_r"])):
        results["results_r"][chr] = np.log2(results["results_r"][chr])

    results["results_r"] = [x.tolist() for x in results["results_r"]]

    for c in range(len(results["results_r"])):
        for i, rR in enumerate(results["results_r"][c]):
            if not np.isfinite(rR):
                results["results_r"][c][i] = 0
                results["results_z"][c][i] = 0
                results["results_w"][c][i] = 0
            if results["results_r"][c][i] != 0:
                results["results_r"][c][i] = (
                    results["results_r"][c][i] - log_r_median
                )


"""
Applies additional blacklist to results_r, results_z
and results_w if requested.
"""


def apply_blacklist(rem_input, results):
    blacklist = _import_bed(rem_input)

    for chr in blacklist.keys():
        for s_e in blacklist[chr]:
            for pos in range(s_e[0], s_e[1]):
                if len(results["results_r"]) < 24 and chr == 23:
                    continue
                if pos >= len(results["results_r"][chr]) or pos < 0:
                    continue
                results["results_r"][chr][pos] = 0
                results["results_z"][chr][pos] = 0
                results["results_w"][chr][pos] = 0


def _import_bed(rem_input):
    bed = {}
    for line in open(rem_input["args"].blacklist):
        chr_name, s, e = line.strip().split("\t")
        if chr_name[:3] == "chr":
            chr_name = chr_name[3:]
        if chr_name == "X":
            chr_name = "23"
        if chr_name == "Y":
            chr_name = "24"
        chr = int(chr_name) - 1
        if chr not in bed.keys():
            bed[chr] = []
        bed[chr].append(
            [
                int(int(s) / rem_input["binsize"]),
                int(int(e) / rem_input["binsize"]) + 1,
            ]
        )
    return bed


def exec_cbs(cna_df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """
    Executed CBS on ratios per bins using weights per bins as weights.

    :param pd.DataFrame cna_df: DataFrame with columns 'chrom', 'maploc', 'null_ratios_mean', 'null_ratios_sd', 'ratios', 'weights' and 'zscores' for each genomic bin
    :param argparse.Namespace args: Command-line arguments namespace
    :return: DataFrame with columns 'sample', 'chrom', 'loc.start', 'loc.end', 'num_mark' and 'seg.mean' for each genomic segment
    :rtype: pd.DataFrame
    """

    # script params
    script_path = Path(__file__).parent / "include" / "CBS.R"

    # Preprocess CNA dataframe for CBS input
    # replace 0s with NaN
    cna_df["ratios"][cna_df["ratios"] == 0] = np.nan
    # replace 0s with a very small number to avoid division by zero in CBS
    cna_df["weights"][cna_df["weights"] == 0] = np.nextafter(0, 1)
    # filter out chromosomes with no valid bins to avoid CBS errors
    cna_df = cna_df.groupby("chrom").filter(
        lambda x: x["ratios"].notna().any()
    )

    # Run CBS and process results
    with TemporaryDirectory() as tmpdir:
        cbs_input = Path(tmpdir) / "cbs_input.json"
        cbs_output = Path(tmpdir) / "cbs_output.json"

        with cbs_input.open("w") as f:
            cna_df.to_json(f, orient="records", indent=4)

        r_cmd = [
            "Rscript",
            str(script_path),
            # IO
            "--infile",
            str(cbs_input),
            "--outfile",
            str(cbs_output),
            # CNA
            "--id",
            str(args.sampleid),
            "--seed",
            str(args.cbs_seed),
            # CBS
            "--alpha",
            str(args.cbs_alpha),
            "--nperm",
            str(args.cbs_nperm),
            "--p_method",
            str(args.cbs_p_method),
            "--min_width",
            str(args.cbs_min_width),
            "--kmax",
            str(args.cbs_kmax),
            "--nmin",
            str(args.cbs_nmin),
            "--eta",
            str(args.cbs_eta),
            "--trim",
            str(args.cbs_trim),
            # set verbose level based on Python logging level (debug -> 3, else 0)
            "--verbose",
            "3" if logging.getLogger().getEffectiveLevel() == 10 else "0",
        ]
        logging.debug("CBS cmd: {}".format(r_cmd))

        try:
            subprocess.check_call(r_cmd)
        except subprocess.CalledProcessError as e:
            logging.critical(f"Rscript failed: {e}")
            cbs_input.rename(os.getcwd() + "/cbs_input.json")
            sys.exit(1)

        with cbs_output.open() as f:
            return pd.read_json(f, orient="records")


def get_segment_z_score(segments_df: pd.DataFrame, cna_df: pd.DataFrame):
    """
    Calculate Z-scores for genomic segments against a baseline.

    For each segment in segment:
        1. Slices background data (ratios, filters, weights) to the segment coordinates.
        2. Filters out bad genomic bins where the filter matrix (results_r) equals 0.
        3. Cleans data by replacing infinite values with NaN.
        4. Computes a column-wise weighted average across the segment slice.
        5. Establishes a null model by calculating the mean and SD of those averages.
        6. Computes the segment's Z-score: (observed_value - null_mean) / null_sd.
        7. Returns "nan" if the baseline mean or SD cannot be mathematically resolved.

    :param pd.DataFrame segments_df: DataFrame with columns 'sample', 'chrom', 'loc.start', 'loc.end', 'num_mark' and 'seg.mean' for each genomic segment
    :param pd.DataFrame cna_df: DataFrame with columns 'chrom', 'maploc', 'null_ratios_mean', 'null_ratios_sd', 'ratios', 'weights' and 'zscores' for each genomic bin

    Returns:
        segments dataframe with an extra column of Z-scores corresponding to each segment in the input dataFrame with invalid calculations marked as "nan".
    """

    # --- STEP 1: Vectorized Slicing via Conditional Merge ---
    # We join segments and cna_filtered on chromosomes, then filter rows
    # where the maploc falls within the segment's start and end boundaries.
    merged = pd.merge(
        segments_df.reset_index().rename(columns={"index": "seg_idx"}),
        cna_df,
        on="chrom",
        suffixes=("_seg", "_cna"),
    )

    # Keep only the genomic bins that physically sit inside each segment's window
    in_window = merged["maploc"].between(
        merged["loc.start"], merged["loc.end"]
    )
    matched_bins = merged[in_window].copy()

    # If no bins matched any segments, append 'nan' column and exit early
    if matched_bins.empty:
        result_df = segments_df.copy()
        result_df["z_score"] = "nan"
        return result_df

    # --- STEP 4 & 5: Vectorized Null Model Math via Groupby ---
    # Pre-calculate weighted products for the weighted mean formula
    matched_bins["weighted_ratios"] = (
        matched_bins["ratios"] * matched_bins["weights"]
    )

    # Group by the unique segment index to calculate mean and SD simultaneously in C-speed
    grouped = matched_bins.groupby("seg_idx")

    stats = grouped.agg(
        total_weights=("weights", "sum"),
        sum_weighted_ratios=("weighted_ratios", "sum"),
        null_mean=("ratios", "mean"),
        null_sd=(
            "ratios",
            lambda x: x.std(ddof=0),
        ),  # Population standard deviation
    )

    # Calculate the localized weighted null mean
    stats["null_weighted_mean"] = np.where(
        stats["total_weights"] > 0,
        stats["sum_weighted_ratios"] / stats["total_weights"],
        np.nan,
    )

    # --- STEP 6, 7 & 8: Compute and Clamp Z-scores Matrix-wide ---
    # Join the computed stats back onto the original segments frame
    result_df = segments_df.copy()
    result_df = result_df.join(stats)

    # Evaluate mathematical validity globally
    invalid_mask = (
        result_df["null_mean"].isna()
        | result_df["null_sd"].isna()
        | (result_df["null_sd"] == 0)
        | result_df["total_weights"].isna()
    )

    # Compute raw Z-scores matrix-wide
    z_raw = (
        result_df["seg.mean"] - result_df["null_weighted_mean"]
    ) / result_df["null_sd"]

    # Apply capping and fill invalid boundaries
    result_df["z_score"] = np.clip(z_raw, -1000, 1000)
    result_df.loc[invalid_mask, "z_score"] = "nan"

    # Drop the temporary calculation helper columns before returning
    columns_to_drop = [
        "total_weights",
        "sum_weighted_ratios",
        "null_mean",
        "null_sd",
        "null_weighted_mean",
    ]
    result_df = result_df.drop(columns=columns_to_drop, errors="ignore")

    return result_df
