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

    :param pd.DataFrame cna_df: DataFrame with columns 'chrom', 'maploc', 'start', 'end', 'null_ratios_mean', 'null_ratios_sd', 'ratios', 'weights' and 'zscores' for each genomic bin
        0 value ratios need to be replaced by NaN
        0 value weights need to be replaced by a very small number to avoid division by zero in CBS
        Assumes dataframe is sorted by chrom and maploc
    :param argparse.Namespace args: Command-line arguments namespace
    :return: DataFrame with columns 'chrom', 'start', 'end', 'num_mark' and 'ratio' for each genomic segment
    :rtype: pd.DataFrame
    """

    # script params
    script_path = Path(__file__).parent / "include" / "CBS.R"

    # Preprocess CNA dataframe for CBS input
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
            segments_df = pd.read_json(f, orient="records")
            segments_df = segments_df.rename(
                columns={
                    "chrom": "chrom",
                    "loc.start": "start",
                    "loc.end": "end",
                    "seg.mean": "ratio",
                }
            )
            return segments_df[["chrom", "start", "end", "ratio"]]


def get_segment_zscore(segments_df: pd.DataFrame, cna_df: pd.DataFrame):
    """
    Calculate Z-scores for genomic segments against a baseline.

    For each segment in segment:
        1. Slices background data (ratios, filters, weights) to the segment coordinates.
        2. Filters out bad genomic bins where the filter matrix (ratios) equals 0.
        3. Cleans data by replacing infinite values with NaN.
        4. Computes a column-wise weighted average across the segment slice.
        5. Establishes a null model by calculating the mean and SD of those averages.
        6. Computes the segment's Z-score: (observed_value - null_mean) / null_sd.
        7. Returns "nan" if the baseline mean or SD cannot be mathematically resolved.

    :param pd.DataFrame segments_df: DataFrame with columns 'chrom', 'start', 'end' and 'ratio' for each genomic segment
    :param pd.DataFrame cna_df: DataFrame with columns 'chrom', 'maploc', 'start', 'end', 'null_ratios_mean', 'null_ratios_sd', 'ratios', 'weights' and 'zscores' for each genomic bin

    Returns:
        segments dataframe with 'chrom', 'start', 'end', 'ratio' and 'zscore' columns with invalid calculations marked as "NaN".
    """
    # Create copies to avoid SettingWithCopyWarning or modifying user inputs unexpectedly
    segments_df_copy = segments_df.copy()
    cna_df_copy = cna_df.copy()

    segments_reset = segments_df_copy.reset_index().rename(
        columns={"index": "seg_idx"}
    )
    merged = pd.merge(
        segments_reset,
        cna_df_copy,
        on="chrom",
        suffixes=("_seg", "_cna"),
    )

    # Keep only the genomic bins that physically sit inside each segment's window
    in_window = merged["maploc"].between(merged["start"], merged["end"])
    matched_bins = merged[in_window].copy()

    # --- STEP 2: Filter out bad genomic bins where the filter matrix (ratios) equals 0 ---
    matched_bins = matched_bins[matched_bins["ratios"] != 0].copy()

    # --- STEP 3: Clean data by replacing infinite values with NaN ---
    cols_to_clean = ["ratios", "weights", "null_ratios_mean", "null_ratios_sd"]
    for col in cols_to_clean:
        if col in matched_bins.columns:
            matched_bins[col] = matched_bins[col].replace(
                [np.inf, -np.inf], np.nan
            )

    if matched_bins.empty:
        segments_df["zscore"] = "nan"
        return segments_df

    # Under the independent bins assumption, combining the null mean and SD for each bin:
    # null_mean = sum(w_i * null_ratios_mean_i) / sum(w_i)
    # null_sd = sqrt(sum(w_i^2 * null_ratios_sd_i^2)) / sum(w_i)
    valid_null = (
        matched_bins["null_ratios_mean"].notna()
        & matched_bins["null_ratios_sd"].notna()
        & matched_bins["weights"].notna()
        & (matched_bins["weights"] > 0)
    )

    matched_bins["w_mean"] = np.where(
        valid_null,
        matched_bins["null_ratios_mean"] * matched_bins["weights"],
        np.nan,
    )
    matched_bins["w_var"] = np.where(
        valid_null,
        (matched_bins["null_ratios_sd"] ** 2) * (matched_bins["weights"] ** 2),
        np.nan,
    )
    matched_bins["w_denom"] = np.where(
        valid_null, matched_bins["weights"], np.nan
    )

    grouped = matched_bins.groupby("seg_idx")
    stats = grouped.agg(
        sum_w_mean=("w_mean", "sum"),
        sum_w_var=("w_var", "sum"),
        sum_w_denom=("w_denom", "sum"),
    )

    stats["null_mean"] = np.where(
        stats["sum_w_denom"] > 0,
        stats["sum_w_mean"] / stats["sum_w_denom"],
        np.nan,
    )
    stats["null_sd"] = np.where(
        stats["sum_w_denom"] > 0,
        np.sqrt(stats["sum_w_var"]) / stats["sum_w_denom"],
        np.nan,
    )

    # Join computed stats back onto the original segments
    result_df = segments_df.copy()
    result_df = result_df.join(stats[["null_mean", "null_sd"]])

    invalid_mask = (
        result_df["null_mean"].isna()
        | result_df["null_sd"].isna()
        | (result_df["null_sd"] == 0)
    )

    # Compute raw Z-scores
    z_raw = (result_df["ratio"] - result_df["null_mean"]) / result_df[
        "null_sd"
    ]

    # Fill invalid boundaries with "nan"
    z_score_series = z_raw.astype(object)
    z_score_series[invalid_mask] = "nan"

    # Assign to segment DataFrames
    segments_df["zscore"] = z_score_series

    return segments_df
