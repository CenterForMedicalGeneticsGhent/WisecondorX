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
import typer
from wisecondorx.predict_control import get_post_processed_result
from typing import Annotated
from wisecondorx.predict_control import normalize
from wisecondorx.predict_output import generate_output_tables
from wisecondorx.plotter import write_plots
from wisecondorx.overall_tools import gender_correct, scale_sample


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


def wcx_predict(
    infile: Path = typer.Argument(..., help=".npz input file"),
    reference: Path = typer.Argument(
        ..., help="Reference .npz, as previously created with newref"
    ),
    outid: str = typer.Argument(
        ...,
        help="Basename (w/o extension) of output files (paths are allowed, e.g. path/to/ID_1)",
    ),
    minrefbins: int = typer.Option(
        150,
        "--minrefbins",
        help="Minimum amount of sensible reference bins per target bin.",
    ),
    maskrepeats: int = typer.Option(
        5,
        "--maskrepeats",
        help="Regions with distances > mean + sd * 3 will be masked. Number of masking cycles.",
    ),
    cbs_alpha: Annotated[
        float,
        typer.Option(
            min=0.0,
            max=1.0,
            help="p-value cut-off for calling a CBS breakpoint.",
        ),
    ] = 1e-4,
    cbs_nperm: Annotated[
        int,
        typer.Option(
            min=0,
            help="Number of permutations for CBS p-value estimation. Higher is more accurate, but also slower.",
        ),
    ] = 10000,
    cbs_p_method: Annotated[
        str,
        typer.Option(
            "--cbs_p_method",
            help="Method for CBS p-value estimation. See DNAcopy documentation for more details.",
        ),
    ] = "hybrid",
    cbs_min_width: Annotated[
        int,
        typer.Option(
            "--cbs_min_width",
            min=2,
            help="Minimum number of bins in a segment for it to be called by CBS. Higher is more stringent.",
        ),
    ] = 2,
    cbs_kmax: Annotated[
        int,
        typer.Option(
            "--cbs_kmax",
            min=1,
            help="Maximum width of smaller segment for permutation in the hybrid method.",
        ),
    ] = 25,
    cbs_nmin: Annotated[
        int,
        typer.Option(
            "--cbs_nmin",
            min=1,
            help="the minimum length of data for which the approximation of maximum statistic is used under the hybrid method. should be larger than 4*kmax",
        ),
    ] = 200,
    cbs_eta: Annotated[
        float,
        typer.Option(
            "--cbs_eta",
            help="the probability to declare a change conditioned on the permuted statistic exceeding the observed statistic exactly j (= 1,...,nperm*alpha) times.",
        ),
    ] = 0.05,
    cbs_trim: Annotated[
        float,
        typer.Option(
            "--cbs_trim",
            help="proportion of data to be trimmed for variance calculation for smoothing outliers and undoing splits based on SD.",
        ),
    ] = 0.025,
    zscore: Annotated[
        float,
        typer.Option(help="z-score cut-off for aberration calling.", min=0.0),
    ] = 5.0,
    beta: Annotated[
        float,
        typer.Option(
            min=0.0,
            max=1.0,
            help="When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations.",
        ),
    ] = None,
    blacklist: Path = typer.Option(
        None,
        "--blacklist",
        help="Blacklist that masks regions in output, structure of header-less file: chr...(/t)startpos(/t)endpos(/n)",
    ),
    gender: Sex = typer.Option(
        None,
        "--gender",
        help="Force WisecondorX to analyze this case as a male (M) or a female (F)",
    ),
    ylim: str = typer.Option(
        "def", "--ylim", help="y-axis limits for plotting. e.g. [-2,2]"
    ),
    bed: bool = typer.Option(
        True,
        "--bed",
        help="Outputs tab-delimited .bed files, containing the most important information",
    ),
    plot: bool = typer.Option(False, "--plot", help="Output .png plots"),
    cairo: bool = typer.Option(
        True,
        "--cairo",
        help="Use cairo for plotting",
    ),
    add_plot_title: bool = typer.Option(
        False,
        "--add-plot-title",
        help="Add the output name as plot title",
    ),
    cbs_seed: int = typer.Option(
        42, "--seed", help="Seed for segmentation algorithm"
    ),
    regions: Path = typer.Option(
        None,
        "--regions",
        help="bed file with regions to be marked on the output plot",
    ),
) -> None:
    """
    Find copy number aberrations.
    """
    logging.info("Starting CNA prediction")

    args = argparse.Namespace()
    args.infile = infile
    args.reference = reference
    args.outid = outid
    args.sampleid = str(os.path.basename(outid))
    args.minrefbins = minrefbins
    args.maskrepeats = maskrepeats
    args.cbs_seed = cbs_seed
    args.cbs_alpha = cbs_alpha
    args.cbs_nperm = cbs_nperm
    args.cbs_p_method = cbs_p_method
    args.cbs_min_width = cbs_min_width
    args.cbs_kmax = cbs_kmax
    args.cbs_nmin = cbs_nmin
    args.cbs_eta = cbs_eta
    args.cbs_trim = cbs_trim
    args.zscore = zscore
    args.beta = beta
    args.blacklist = blacklist
    args.gender = gender
    args.ylim = ylim
    args.bed = bed
    args.plot = plot
    args.cairo = cairo
    args.add_plot_title = add_plot_title
    args.regions = regions

    if not args.bed and not args.plot:
        logging.critical(
            "No output format selected. "
            "Select at least one of the supported output formats (--bed, --plot)"
        )
        sys.exit()

    if args.zscore <= 0:
        logging.critical(
            "Parameter --zscore should be a strictly positive number"
        )
        sys.exit()

    if args.beta is not None:
        if args.beta <= 0 or args.beta > 1:
            logging.critical(
                "Parameter --beta should be a strictly positive number lower than or equal to 1"
            )
            sys.exit()

    if args.cbs_alpha <= 0 or args.cbs_alpha > 1:
        logging.critical(
            "Parameter --cbs_alpha should be a strictly positive number lower than or equal to 1"
        )
        sys.exit()

    logging.info("Importing data ...")
    ref_file = np.load(args.reference, encoding="latin1", allow_pickle=True)
    sample_file = np.load(args.infile, encoding="latin1", allow_pickle=True)

    sample = sample_file["sample"].item()
    n_reads = sum([sum(sample[x]) for x in sample.keys()])

    sample = scale_sample(
        sample, int(sample_file["binsize"].item()), int(ref_file["binsize"])
    )

    gender = (
        args.gender
        if args.gender
        else predict_gender(sample, ref_file["trained_cutoff"])
    )
    if not ref_file["is_nipt"]:
        sample = gender_correct(sample, gender)
        ref_gender = gender
    else:
        ref_gender = Sex.FEMALE

    logging.info("Normalizing autosomes ...")

    results_r, results_z, results_w, ref_sizes, m_lr, m_z = normalize(
        args, sample, ref_file, "A"
    )

    if not ref_file["is_nipt"]:
        if not ref_file["has_male"] and gender == Sex.MALE:
            logging.warning(
                "This sample is male, whilst the reference is created with fewer than 5 males. "
                "The female gonosomal reference will be used for X predictions. Note that these might "
                "not be accurate. If the latter is desired, create a new reference and include more "
                "male samples."
            )
            ref_gender = Sex.FEMALE

        elif not ref_file["has_female"] and gender == Sex.FEMALE:
            logging.warning(
                "This sample is female, whilst the reference is created with fewer than 5 females. "
                "The male gonosomal reference will be used for XY predictions. Note that these might "
                "not be accurate. If the latter is desired, create a new reference and include more "
                "female samples."
            )
            ref_gender = Sex.MALE

    logging.info("Normalizing gonosomes ...")

    null_ratios_aut_per_bin = ref_file["null_ratios"]
    null_ratios_gon_per_bin = ref_file[
        "null_ratios.{}".format(ref_gender.value)
    ][len(null_ratios_aut_per_bin) :]

    results_r_2, results_z_2, results_w_2, ref_sizes_2, _, _ = normalize(
        args, sample, ref_file, ref_gender.value
    )

    rem_input = {
        "args": args,
        "wd": str(os.path.dirname(os.path.realpath(__file__))),
        "binsize": int(ref_file["binsize"]),
        "n_reads": n_reads,
        "ref_gender": ref_gender.value,
        "gender": gender.value,
        "mask": ref_file["mask.{}".format(ref_gender.value)],
        "bins_per_chr": ref_file["bins_per_chr.{}".format(ref_gender.value)],
        "masked_bins_per_chr": ref_file[
            "masked_bins_per_chr.{}".format(ref_gender.value)
        ],
        "masked_bins_per_chr_cum": ref_file[
            "masked_bins_per_chr_cum.{}".format(ref_gender.value)
        ],
    }

    del ref_file

    results_r = np.append(results_r, results_r_2)
    results_z = np.append(results_z, results_z_2) - m_z
    results_w = np.append(
        results_w * np.nanmean(results_w_2),
        results_w_2 * np.nanmean(results_w),
    )
    results_w = results_w / np.nanmean(results_w)

    if np.isnan(results_w).any() or np.isinf(results_w).any():
        logging.warning(
            "Non-numeric values found in weights -- reference too small. Circular binary segmentation and z-scoring will be unweighted"
        )
        results_w = np.ones(len(results_w))

    ref_sizes = np.append(ref_sizes, ref_sizes_2)

    null_ratios = np.array(
        [x.tolist() for x in null_ratios_aut_per_bin]
        + [x.tolist() for x in null_ratios_gon_per_bin],
        dtype=object,
    )

    results = {
        "results_r": results_r,
        "results_z": results_z,
        "results_w": results_w,
        "results_nr": null_ratios,
    }

    for result in results.keys():
        results[result] = get_post_processed_result(
            args, results[result], ref_sizes, rem_input
        )

    log_trans(results, m_lr)

    if args.blacklist:
        logging.info("Applying blacklist ...")
        apply_blacklist(rem_input, results)

    ## Vectorize results for easier handling in post-processing and plotting
    # vector containing chromosome number for each bin throughout the genome
    chrom_vec = np.repeat(
        list(range(1, len(results["results_r"]) + 1)),
        [len(c) for c in results["results_r"]],
    )
    # vector containing bin number for each bin throughout the genome
    maploc_vec = np.concatenate(
        [np.arange(1, n + 1) for n in [len(c) for c in results["results_r"]]]
    )
    # vector containing start position for each bin throughout the genome
    start_pos_vec = (maploc_vec - 1) * rem_input["binsize"]
    # vector containing end position for each bin throughout the genome
    end_pos_vec = maploc_vec * rem_input["binsize"]
    # vector containing ratio for each bin
    ratios_vec = np.concatenate(results["results_r"]).astype(float)
    # replace 0s with NaN
    ratios_vec[ratios_vec == 0] = np.nan

    # vector containing z-score for each bin
    zscores_vec = np.concatenate(results["results_z"]).astype(float)

    # vector containing weight for each bin
    weights_vec = np.concatenate(results["results_w"]).astype(float)
    # replace 0s with a very small number to avoid division by zero in CBS
    weights_vec[weights_vec == 0] = np.nextafter(0, 1)
    if np.isnan(weights_vec).any() or np.isinf(weights_vec).any():
        logging.warning(
            "Non-numeric values found in weights -- reference too small. Circular binary segmentation and z-scoring will be unweighted"
        )
        weights_vec = np.ones(len(weights_vec))

    # vectors containing mean and stdev of null ratios for each bin
    null_ratios_mean_vec = np.array(
        [
            float(np.mean(bin_null_ratios))
            for chr_null_ratios in results["results_nr"]
            for bin_null_ratios in chr_null_ratios
        ],
        dtype=float,
    )
    null_ratios_stdev_vec = np.array(
        [
            float(np.std(bin_null_ratios))
            for chr_null_ratios in results["results_nr"]
            for bin_null_ratios in chr_null_ratios
        ],
        dtype=float,
    )

    cna_df = (
        pd.DataFrame(
            {
                "chrom": chrom_vec,
                "maploc": maploc_vec,
                "start": start_pos_vec,
                "end": end_pos_vec,
                "null_ratios_mean": null_ratios_mean_vec,
                "null_ratios_sd": null_ratios_stdev_vec,
                "ratios": ratios_vec,
                "weights": weights_vec,
                "zscores": zscores_vec,
            }
        )
        .sort_values(by=["chrom", "maploc"], ascending=[True, True])
        .reset_index(drop=True)
    )

    logging.info("Executing circular binary segmentation ...")

    # Run CBS and add z-scores per segment
    segments_df = get_segment_zscore(exec_cbs(cna_df, args), cna_df)
    # Rename columns, replace start and end idx with genomic positions, and sort by chromosome and start position
    segments_df["start"] = segments_df["start"] * rem_input["binsize"]
    segments_df["end"] = segments_df["end"] * rem_input["binsize"]
    segments_df = segments_df.sort_values(
        by=["chrom", "start"], ascending=[True, True]
    ).reset_index(drop=True)

    if args.bed:
        logging.info("Writing tables ...")
        generate_output_tables(cna_df, segments_df, args)

    if args.plot:
        logging.info("Writing plots ...")
        data = {
            "ref_gender": str(rem_input["ref_gender"]),
            "beta": str(args.beta),
            "zscore": str(args.zscore),
            "binsize": str(rem_input["binsize"]),
            "n_reads": str(rem_input["n_reads"]),
            "cairo": args.cairo,
            "cna_df": cna_df,
            "segments_df": segments_df,
            "ylim": str(args.ylim),
            "regions": str(args.regions),
            "out_dir": str(f"{args.outid}.plots"),
        }

        if args.add_plot_title:
            data["plot_title"] = str(args.sampleid)

        write_plots(data)

    logging.info("Finished prediction")


def wcx_gender(
    infile: Path = typer.Argument(..., help=".npz input file"),
    reference: Path = typer.Argument(
        ..., help="Reference .npz, as previously created with newref"
    ),
) -> None:
    """
    Returns the sex of a .npz resulting from convert, based on a Gaussian mixture model trained during the newref phase.
    """

    ref_file = np.load(reference, encoding="latin1", allow_pickle=True)
    sample_file = np.load(infile, encoding="latin1", allow_pickle=True)
    gender = predict_gender(
        sample_file["sample"].item(), ref_file["trained_cutoff"]
    )
    if gender == "M":
        print("male")
    else:
        print("female")
