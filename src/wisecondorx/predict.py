from typing import Any
import json
import subprocess
import math
import sys
from pathlib import Path
from typing import Annotated
import numpy as np
import pandas as pd
from scipy.stats import norm
import typer
from dataclasses import dataclass
import re
from tempfile import TemporaryDirectory
from wisecondorx.utils import (
    scale_bins_per_chromosome,
    sex_correct,
    Sex,
)
import os
from wisecondorx.plotter import write_plots

import logging


def wcx_sex(
    infile: Path = typer.Argument(..., help=".npz input file"),
    reference: Path = typer.Argument(
        ..., help="Reference .npz, as previously created with newref"
    ),
) -> None:
    """
    Returns the sex of a .npz resulting from convert, based on a Gaussian mixture model trained during the newref phase.
    """
    reference_df = np.load(reference, encoding="latin1", allow_pickle=True)
    sample_df = np.load(infile, encoding="latin1", allow_pickle=True)
    sample_sex = predict_sex(
        sample_df["sample"].item(), reference_df["trained_cutoff"]
    )
    if sample_sex == Sex.MALE:
        print("Male")
    elif sample_sex == Sex.FEMALE:
        print("Female")
    else:
        print("Unknown")


def wcx_predict(
    infile: Path = typer.Argument(..., help=".npz input file"),
    reference: Path = typer.Argument(
        ..., help="Reference .npz, as previously created with newref"
    ),
    prefix: Path = typer.Argument(
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
    alpha: Annotated[
        float,
        typer.Option(
            min=0.0,
            max=1.0,
            help="p-value cut-off for calling a CBS breakpoint.",
        ),
    ] = 1e-4,
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
    sex_override: Sex = typer.Option(
        None,
        "--sex",
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
    add_plot_title: bool = typer.Option(
        False,
        "--add-plot-title",
        help="Add the output name as plot title",
    ),
    seed: int = typer.Option(
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

    if not bed and not plot:
        logging.critical(
            "No output format selected. Select at least one of the supported output formats (--bed, --plot)"
        )
        sys.exit(1)

    logging.info("Importing data ...")
    # Unpack reference npz
    # Reference npz keys
    # 'has_female', 'has_male',
    # 'binsize', 'mask',
    # 'bins_per_chr', 'masked_bins_per_chr' (count of masked bins per chromosome), 'masked_bins_per_chr_cum',
    # 'pca_components', 'pca_mean',
    # 'indexes', 'distances', 'null_ratios',
    # 'binsize.F', 'mask.F', 'bins_per_chr.F', 'masked_bins_per_chr.F', 'masked_bins_per_chr_cum.F', 'pca_components.F', 'pca_mean.F', 'indexes.F', 'distances.F', 'null_ratios.F',
    # 'binsize.M', 'mask.M', 'bins_per_chr.M', 'masked_bins_per_chr.M', 'masked_bins_per_chr_cum.M', 'pca_components.M', 'pca_mean.M', 'indexes.M', 'distances.M', 'null_ratios.M',
    # 'is_nipt', 'trained_cutoff'
    reference_npz = np.load(reference, encoding="latin1", allow_pickle=True)
    reference_binsize: int = reference_npz["binsize"]
    reference_is_nipt = reference_npz["is_nipt"]
    reference_trained_cutoff = reference_npz["trained_cutoff"]

    # Unpack sample npz
    sample_npz = np.load(infile, encoding="latin1", allow_pickle=True)
    sample_bins_per_chr: dict[str, np.ndarray] = sample_npz.get(
        "sample", sample_npz.get("reads_per_bin")
    ).item()
    sample_binsize = sample_npz["binsize"]
    sample_total_reads = sum(
        [sum(sample_bins_per_chr[chr]) for chr in sample_bins_per_chr.keys()]
    )
    scaled_sample_bins_per_chromosome = scale_bins_per_chromosome(
        sample_bins_per_chr, sample_binsize, reference_binsize
    )
    sample_sex = (
        sex_override
        if sex_override
        else predict_sex(
            scaled_sample_bins_per_chromosome, reference_trained_cutoff
        )
    )
    if reference_is_nipt:
        reference_sex = Sex.FEMALE
    else:
        scaled_sample_bins_per_chromosome = sex_correct(
            scaled_sample_bins_per_chromosome, sample_sex
        )
        reference_sex = sample_sex
        if not reference_npz["has_male"] and sample_sex == Sex.MALE:
            logging.warning(
                "This sample is male, whilst the reference is created with fewer than 5 males. The female gonosomal reference will be used for X predictions."
            )
            reference_sex = Sex.FEMALE
        elif not reference_npz["has_female"] and sample_sex == Sex.FEMALE:
            logging.warning(
                "This sample is female, whilst the reference is created with fewer than 5 females. The male gonosomal reference will be used for XY predictions."
            )
            reference_sex = Sex.MALE

    logging.info("Normalizing autosomes ...")
    (
        autosomal_ratios,
        autosomal_zscores,
        autosomal_weights,
        autosomal_ref_sizes,
        median_logratio,
        median_zscore,
    ) = normalize(
        sample=scaled_sample_bins_per_chromosome,
        reference_npz=reference_npz,
        ref_sex=Sex.AUTOSOMAL,
        maskrepeats=maskrepeats,
    )

    logging.info("Normalizing gonosomes ...")
    (
        gonosomal_ratios,
        gonosomal_zscores,
        gonosomal_weights,
        gonosomal_ref_sizes,
        _,
        _,
    ) = normalize(
        sample=scaled_sample_bins_per_chromosome,
        reference_npz=reference_npz,
        ref_sex=reference_sex,
        maskrepeats=maskrepeats,
    )

    rem_input = {
        "binsize": reference_binsize,
        "n_reads": sample_total_reads,
        "ref_sex": reference_sex,
        "sex": reference_sex,
        "mask": reference_npz[f"mask.{reference_sex.value}"],
        "bins_per_chr": reference_npz[f"bins_per_chr.{reference_sex.value}"],
        "masked_bins_per_chr": reference_npz[
            f"masked_bins_per_chr.{reference_sex.value}"
        ],
        "masked_bins_per_chr_cum": reference_npz[
            f"masked_bins_per_chr_cum.{reference_sex.value}"
        ],
    }

    ref_sizes = np.append(autosomal_ref_sizes, gonosomal_ref_sizes)
    ratios = np.append(autosomal_ratios, gonosomal_ratios)
    zscores = np.append(autosomal_zscores, gonosomal_zscores) - median_zscore
    weights = np.append(
        autosomal_weights * np.nanmean(gonosomal_weights),
        gonosomal_weights * np.nanmean(autosomal_weights),
    )
    weights = weights / np.nanmean(weights)

    if np.isnan(weights).any() or np.isinf(weights).any():
        logging.warning(
            "Non-numeric values found in weights -- reference too small. Circular binary segmentation and z-scoring will be unweighted"
        )
        weights = np.ones(len(weights))

    null_ratios_aut_per_bin = reference_npz["null_ratios"]
    null_ratios_gon_per_bin = reference_npz[
        f"null_ratios.{reference_sex.value}"
    ][len(null_ratios_aut_per_bin) :]

    null_ratios = np.array(
        [x.tolist() for x in null_ratios_aut_per_bin]
        + [x.tolist() for x in null_ratios_gon_per_bin],
        dtype=object,
    )

    results = {
        "ratios": ratios,
        "zscores": zscores,
        "weights": weights,
        "null_ratios": null_ratios,
    }

    for key, result in results.items():
        results[key] = get_post_processed_result(
            minrefbins=minrefbins,
            result=result,
            ref_sizes=ref_sizes,
            rem_input=rem_input,
        )

    logtransform(results, median_logratio)

    if blacklist:
        logging.info("Applying blacklist ...")
        results = apply_blacklist(
            blacklist_path=blacklist,
            binsize=reference_binsize,
            results=results,
        )

    logging.info("Executing circular binary segmentation ...")

    # segments is a list of lists, where each inner list is a segment
    # each segment is a list of 4 values: [chr, start, end, log2_ratio]
    segments = run_cbs(
        alpha=alpha,
        seed=seed,
        ref_sex=reference_sex,
        binsize=reference_binsize,
        ratios=results["ratios"],
        weights=results["weights"],
    )

    segment_zscores = get_z_score(
        segments,
        results["ratios"],
        results["null_ratios"],
        results["weights"],
    )

    for i in range(len(segments)):
        segments[i].append(segment_zscores[i])
    results["segments"] = segments
    if bed:
        logging.info("Writing tables ...")
        write_tables(
            prefix=prefix,
            binsize=reference_binsize,
            regions=regions,
            bins_per_chr=scaled_sample_bins_per_chromosome,
            ref_sex=reference_sex,
            beta=beta,
            zscore=zscore,
            sex=sample_sex,
            n_reads=sample_total_reads,
            ratios=ratios,
            weights=weights,
            segments=segments,
        )

    if plot:
        logging.info("Writing plots ...")
        write_plots(
            out_dir=Path(f"{prefix}.plots"),
            ref_sex=reference_sex,
            beta=beta,
            zscore=zscore,
            binsize=reference_binsize,
            n_reads=sample_total_reads,
            ylim_str=ylim,
            regions_file=regions,
            plot_title=str(prefix.name) if add_plot_title else "",
            ratios=ratios,
            weights=weights,
            segments=segments,
        )

    logging.info("Finished prediction")


def predict_sex(sample: dict[str, np.ndarray], trained_cutoff: float) -> Sex:
    """
    Returns sex based on Gaussian mixture
    model trained during newref phase.
    """
    Y_fraction = float(np.sum(sample["24"])) / float(
        np.sum([np.sum(sample[x]) for x in sample.keys()])
    )
    if Y_fraction > trained_cutoff:
        return Sex.MALE
    else:
        return Sex.FEMALE


def inflate_results(results, rem_input) -> list[float]:
    """
    Unmasks results array.
    """
    temp = [0 for x in rem_input["mask"]]
    j = 0
    for i, val in enumerate(rem_input["mask"]):
        if val:
            temp[i] = results[j]
            j += 1
    return temp


def logtransform(
    results: dict[str, list[list[float]]], log_r_median: float
) -> dict[str, list[list[float]]]:
    """
    Log2-transforms ratios. If resulting elements are infinite,
    all corresponding possible positions (at ratios, zscores
    and weights are set to 0 (blacklist)).
    """
    for chr in range(len(results["ratios"])):
        results["ratios"][chr] = np.log2(results["ratios"][chr])

    results["ratios"] = [x.tolist() for x in results["ratios"]]

    for c in range(len(results["ratios"])):
        for i, rR in enumerate(results["ratios"][c]):
            if not np.isfinite(rR):
                results["ratios"][c][i] = 0
                results["zscores"][c][i] = 0
                results["weights"][c][i] = 0
            if results["ratios"][c][i] != 0:
                results["ratios"][c][i] = (
                    results["ratios"][c][i] - log_r_median
                )
    return results


def apply_blacklist(
    blacklist_path: Path, binsize: int, results: dict[str, Any]
) -> dict[str, list[list[float]]]:
    """
    Applies additional blacklist to ratios, zscores
    and weights if requested.
    """
    blacklist = _import_bed(blacklist_path, binsize)

    for chr in blacklist.keys():
        for s_e in blacklist[chr]:
            for pos in range(s_e[0], s_e[1]):
                if len(results["ratios"]) < 24 and chr == 23:
                    continue
                if pos >= len(results["ratios"][chr]) or pos < 0:
                    continue
                results["ratios"][chr][pos] = 0
                results["zscores"][chr][pos] = 0
                results["weights"][chr][pos] = 0
    return results


def _import_bed(
    blacklist_path: Path, binsize: int
) -> dict[int, list[list[int]]]:
    """
    Imports a blacklist from a BED file.
    """
    bed = {}
    for line in open(blacklist_path):
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
                int(int(s) / binsize),
                int(int(e) / binsize) + 1,
            ]
        )
    return bed


def normalize(sample, reference_npz, ref_sex, maskrepeats):
    if ref_sex == Sex.AUTOSOMAL:
        suffix = ""
        chrom_idx = 0
        chrom_total = 0
    else:
        suffix = ".{}".format(ref_sex.value)
        chrom_idx = 22
        chrom_total = reference_npz[
            "masked_bins_per_chr_cum{}".format(suffix)
        ][chrom_idx - 1]

    sample = coverage_normalize_and_mask(
        sample,
        reference_npz[f"bins_per_chr{suffix}"],
        reference_npz[f"mask{suffix}"],
    )
    sample = project_pc(
        sample,
        reference_npz[f"pca_components{suffix}"],
        reference_npz[f"pca_mean{suffix}"],
    )
    weights = get_weights(reference_npz["distances{}".format(suffix)])[
        chrom_total:
    ]
    optimal_cutoff = get_optimal_cutoff(
        reference_npz["distances"], maskrepeats
    )
    zscores, ratios, ref_sizes, median_logratio, median_zscore = (
        normalize_repeat(
            sample,
            reference_npz,
            optimal_cutoff,
            chrom_total,
            chrom_idx,
            suffix,
        )
    )

    return ratios, zscores, weights, ref_sizes, median_logratio, median_zscore


def coverage_normalize_and_mask(
    sample, reference_bins_per_chr, reference_mask
) -> np.ndarray:
    by_chr = []

    chrs = range(1, len(reference_bins_per_chr) + 1)

    for chr in chrs:
        this_chr = np.zeros(reference_bins_per_chr[chr - 1], dtype=float)
        min_len = min(reference_bins_per_chr[chr - 1], len(sample[str(chr)]))
        this_chr[:min_len] = sample[str(chr)][:min_len]
        by_chr.append(this_chr)

    all_chr = np.concatenate(by_chr, axis=0)
    all_chr = all_chr / np.sum(all_chr)
    masked_data = all_chr[reference_mask]

    return masked_data


def project_pc(sample, pca_components, pca_mean):
    """
    Project test sample to PCA space.
    """
    # Center the data
    centered = sample - pca_mean
    # Project to PCA space (equivalent to pca.transform)
    transform = np.dot(centered, pca_components.T)
    # Reconstruct from PCA space (equivalent to pca.inverse_transform)
    reconstructed = np.dot(transform, pca_components) + pca_mean
    return sample / reconstructed


def get_weights(reference_distances: np.ndarray) -> np.ndarray:
    """
    The means of sets of within-sample reference
    distances can serve as inverse weights for
    CBS, Z-scoring and plotting.
    """
    inverse_weights = [np.mean(np.sqrt(x)) for x in reference_distances]
    weights = np.array([1 / x for x in inverse_weights])
    return weights


def get_optimal_cutoff(reference_distances: np.ndarray, repeats: int) -> float:
    """
    Defines cutoff that will add bins to a blacklist
    depending on the within reference distances.
    """
    cutoff = float("inf")
    for i in range(0, repeats):
        mask = reference_distances < cutoff
        average = np.average(reference_distances[mask])
        stddev = np.std(reference_distances[mask])
        cutoff = average + 3 * stddev
    return cutoff


def normalize_repeat(
    sample, reference_npz, optimal_cutoff, chr_total, chr_idx, suffix
):
    """
    Within sample normalization. Cycles through a number
    of repeats where z-scores define whether a bin is seen
    as 'normal' in a sample (in most cases this means
    'non-aberrant') -- if it is, it can be used as a
    reference for the other bins.
    """
    zscores = None
    ratios = None
    ref_sizes = None
    sample_copy = np.copy(sample)
    for i in range(3):
        zscores, ratios, ref_sizes = _normalize_once(
            sample,
            sample_copy,
            reference_npz,
            optimal_cutoff,
            chr_total,
            chr_idx,
            suffix,
        )

        sample_copy[chr_total:][np.abs(zscores) >= norm.ppf(0.99)] = -1
    median_logratio = np.nanmedian(np.log2(ratios))
    median_zscore = np.nanmedian(zscores)

    return zscores, ratios, ref_sizes, median_logratio, median_zscore


def _normalize_once(
    sample, sample_copy, ref_file, optimal_cutoff, chr_total, chr_idx, suffix
):
    masked_bins_per_chr = ref_file["masked_bins_per_chr{}".format(suffix)]
    masked_bins_per_chr_cum = ref_file[
        "masked_bins_per_chr_cum{}".format(suffix)
    ]
    zscores = np.zeros(masked_bins_per_chr_cum[-1])[chr_total:]
    ratios = np.zeros(masked_bins_per_chr_cum[-1])[chr_total:]
    ref_sizes = np.zeros(masked_bins_per_chr_cum[-1])[chr_total:]
    indexes = ref_file["indexes{}".format(suffix)]
    distances = ref_file["distances{}".format(suffix)]

    i = chr_total
    i2 = 0
    for chr in list(range(len(masked_bins_per_chr)))[chr_idx:]:
        start = masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]
        end = masked_bins_per_chr_cum[chr]
        chr_data = np.concatenate(
            (
                sample_copy[
                    : masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]
                ],
                sample_copy[masked_bins_per_chr_cum[chr] :],
            )
        )

        for index in indexes[start:end]:
            ref_data = chr_data[index[distances[i] < optimal_cutoff]]
            ref_data = ref_data[ref_data >= 0]
            ref_stdev = np.std(ref_data)
            zscores[i2] = (sample[i] - np.mean(ref_data)) / ref_stdev
            ratios[i2] = sample[i] / np.median(ref_data)
            ref_sizes[i2] = ref_data.shape[0]
            i += 1
            i2 += 1

    return zscores, ratios, ref_sizes


def get_post_processed_result(
    minrefbins: int,
    ref_sizes,
    rem_input,
    result,
):
    """
    Function processes a result (e.g. ratios)
    to an easy-to-interpret format. Bins without
    information are set to 0.
    """

    infinite_mask = ref_sizes < minrefbins
    result[infinite_mask] = 0
    inflated_results = inflate_results(result, rem_input)

    final_results = []
    for chr in range(len(rem_input["bins_per_chr"])):
        chr_data = inflated_results[
            sum(rem_input["bins_per_chr"][:chr]) : sum(
                rem_input["bins_per_chr"][: chr + 1]
            )
        ]
        final_results.append(chr_data)

    return final_results


"""
Calculates zz-scores, marks aberrations and
writes tables.
"""


def write_tables(
    prefix: Path,
    binsize: int,
    regions: str,
    bins_per_chr: list,
    ref_sex: Sex,
    beta: float,
    zscore: float,
    sex: Sex,
    n_reads: int,
    results,
) -> None:
    write_bins_bed(prefix, binsize, results)
    write_segments_and_aberrations_bed(
        prefix, binsize, ref_sex, beta, zscore, results
    )
    write_statistics_file(prefix, bins_per_chr, binsize, sex, n_reads, results)
    if regions is not None:
        write_regions_bed(prefix, binsize, regions, bins_per_chr, results)


def write_bins_bed(prefix: Path, binsize: int, results) -> None:
    bins_file = open(Path(f"{prefix}_bins.bed"), "w")
    bins_file.write("chr\tstart\tend\tid\tratio\tzscore\n")
    ratios = results["ratios"]
    zscores = results["zscores"]
    zscores = results["zscores"]

    for chr in range(len(ratios)):
        chr_name = str(chr + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"
        feat = 1
        for i in range(len(ratios[chr])):
            r = ratios[chr][i]
            z = zscores[chr][i]
            if r == 0:
                r = "nan"
            if z == 0:
                z = "nan"
            feat_str = "{}:{}-{}".format(
                chr_name, str(feat), str(feat + binsize - 1)
            )
            row = [
                chr_name,
                feat,
                feat + binsize - 1,
                feat_str,
                round(float(r), 4),
                round(float(z), 4),
            ]
            bins_file.write("{}\n".format("\t".join([str(x) for x in row])))
            feat += binsize
    bins_file.close()


def write_regions_bed(
    prefix: Path,
    binsize: int,
    regions_path: str,
    bins_per_chr: list[int],
    results,
) -> None:
    regions_file = open(Path(f"{prefix}_regions.bed"), "w")
    regions_file.write("chr\tstart\tend\tname\tratio\tzscore\n")

    with open(regions_path, "r") as regions_file_handle:
        regions = [
            line.strip().split("\t")
            for line in regions_file_handle
            if line.strip() != ""
        ]

        for region in regions:
            assert len(region) >= 4, (
                "Regions file must have at least 4 columns: chr, start, end, name"
            )
            chr_name, start, end, name = (
                region[0],
                region[1],
                region[2],
                region[3],
            )

            # Convert chromosome name to zero-based index
            if chr_name == "chrX" or chr_name == "X":
                chr = 21
            if chr_name == "chrY" or chr_name == "Y":
                chr = 22
            chr = int(re.sub("chr", "", chr_name)) - 1
            start_bin = int(start) // binsize
            end_bin = int(end) // binsize
            if end_bin >= bins_per_chr[chr]:
                end_bin = bins_per_chr[chr] - 1

            if start_bin < 0 or end_bin < 0 or start_bin > end_bin:
                regions_file.write(
                    "Skipping invalid region: {}\n".format("\t".join(region))
                )
                continue

            # Extract ratios, weights, and z-scores for the region
            region_ratios = results["ratios"][chr][start_bin : end_bin + 1]
            region_weights = results["weights"][chr][start_bin : end_bin + 1]
            region_zscores = results["zscores"][chr][start_bin : end_bin + 1]

            if len(region_ratios) == 0:
                regions_file.write(
                    "Skipping region with no bins: {}\n".format(
                        "\t".join(region)
                    )
                )
                continue

            # Calculate weighted means
            ratio_mean = np.ma.average(region_ratios, weights=region_weights)
            zscore_mean = np.ma.average(region_zscores, weights=region_weights)

            if ratio_mean == 0:
                ratio_mean = "nan"
            if zscore_mean == 0:
                zscore_mean = "nan"

            row = [chr_name, start, end, name, ratio_mean, zscore_mean]
            regions_file.write("{}\n".format("\t".join([str(x) for x in row])))

    regions_file.close()


def write_segments_and_aberrations_bed(
    prefix: Path,
    binsize: int,
    ref_sex: Sex,
    beta: float,
    zscore: float,
    results,
) -> None:
    segments_file = open(Path(f"{prefix}_segments.bed"), "w")
    aberrations_file = open(Path(f"{prefix}_aberrations.bed"), "w")
    segments_file.write("chr\tstart\tend\tratio\tzscore\n")
    aberrations_file.write("chr\tstart\tend\tratio\tzscore\ttype\n")

    for segment in results["segments"]:
        chr_name = str(segment[0] + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"
        row = [
            chr_name,
            int(segment[1] * binsize + 1),
            int(segment[2] * binsize),
            segment[4],
            segment[3],
        ]
        segments_file.write("{}\n".format("\t".join([str(x) for x in row])))

        ploidy = 2
        if (chr_name == "X" or chr_name == "Y") and ref_sex == Sex.MALE:
            ploidy = 1

        if beta:
            loss_cutoff = float(np.log2((ploidy - (beta / 2)) / ploidy))
            gain_cutoff = float(np.log2((ploidy + (beta / 2)) / ploidy))
            if float(segment[4]) > gain_cutoff:
                aberrations_file.write(
                    "{}\tgain\n".format("\t".join([str(x) for x in row]))
                )
            elif float(segment[4]) < loss_cutoff:
                aberrations_file.write(
                    "{}\tloss\n".format("\t".join([str(x) for x in row]))
                )
        elif isinstance(segment[3], str):
            continue
        else:
            if float(segment[3]) > zscore:
                aberrations_file.write(
                    "{}\tgain\n".format("\t".join([str(x) for x in row]))
                )
            elif float(segment[3]) < -zscore:
                aberrations_file.write(
                    "{}\tloss\n".format("\t".join([str(x) for x in row]))
                )

    segments_file.close()
    aberrations_file.close()


@dataclass
class ChrStats:
    chr: str
    ratio_mean: float
    ratio_median: float
    zscore: float


def write_statistics_file(
    prefix: Path,
    bins_per_chr: dict[str, np.ndarray],
    binsize: int,
    sex: Sex,
    n_reads: int,
    ratios,
    segments,
    weights,
) -> None:
    chr_stats: list[ChrStats] = []
    for chr in range(len(ratios)):
        chr_name = str(chr + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"

        ratio_mean = np.ma.average(ratios[chr], weights=weights[chr])
        ratio_median = np.median([x for x in ratios[chr] if x != 0])
        # result_c = [chr, 0, bins_per_chr[chr] - 1, ratio_mean]
        # zscore = get_z_score(result_c, results)[0]
        zscore = None

        chr_stats.append(ChrStats(chr_name, ratio_mean, ratio_median, zscore))

    stats_json = open(Path(f"{prefix}_statistics.json"), "w")
    stats_dict = {
        "sex": sex,
        "n_reads": n_reads,
        "chr_stats": chr_stats,
        "ratio_stdev": round(
            float(np.nanstd([c.ratio_mean for c in chr_stats])), 4
        ),
        "msv": round(
            get_median_segment_variance(segments, ratios),
            4,
        ),
        "cpa": round(get_cpa(segments, binsize), 4),
    }
    json.dump(stats_dict, stats_json, indent=4)
    stats_json.close()


def get_cpa(segments, binsize: int) -> float:
    """
    Returns CPA, measure for sample-wise abnormality.
    """
    x: float = 0.0
    for segment in segments:
        x += (segment[2] - segment[1] + 1) * binsize * abs(segment[3])
    cpa = x / len(segments) * (10**-8)
    return cpa


def get_median_segment_variance(segments, ratios) -> float:
    """
    Returns median segment variance, measure for sample-wise noise.
    """
    vars: list[float] = []
    for segment in segments:
        segment_chr, segment_start, segment_stop, segment_ratio = segment
        segment_ratio = ratios[segment_chr][segment_start:segment_stop]
        segment_ratio = [ratio for ratio in segment_ratio if ratio != 0]
        if segment_ratio:
            var = np.var(segment_ratio)
            vars.append(var)
    return np.median(vars)


def run_cbs(
    alpha: float,
    seed: int,
    ref_sex: Sex,
    binsize: int,
    ratios_per_chr: np.ndarray,
    weights_per_chr: np.ndarray,
):
    script_path = Path(__file__).parent / "include" / "CBS.R"
    segment_ratios = []
    with TemporaryDirectory() as tmpdir:
        cbs_input = Path(tmpdir) / "cbs_input.json"
        cbs_output = Path(tmpdir) / "cbs_output.json"

        json.dump(
            {
                "ratios": ratios_per_chr,
                "weights": weights_per_chr,
            },
            open(cbs_input, "w"),
            indent=4,
        )

        r_cmd = [
            "Rscript",
            script_path,
            "--infile",
            cbs_input,
            "--outfile",
            cbs_output,
            "--sex",
            ref_sex.value,
            "--alpha",
            str(alpha),
            "--binsize",
            str(binsize),
            "--seed",
            str(seed),
        ]
        logging.debug(f"CBS cmd: {r_cmd}")

        try:
            subprocess.check_call(r_cmd)
            segment_ratios = _get_processed_cbs(cbs_output)
        except subprocess.CalledProcessError as e:
            logging.critical(f"Rscript failed: {e}")
            cbs_input.rename(os.getcwd() + "/cbs_input.json")
            sys.exit(1)

    return segment_ratios


def _get_processed_cbs(
    cbs_json: Path,
) -> list[tuple[int, int, int, float]]:
    segments = []
    with cbs_json.open("r") as f:
        cbs_data = json.load(f)
    for segment in cbs_data:
        chr = int(segment["chromosome"]) - 1
        start = int(segment["start"])
        end = int(segment["end"])
        ratio = segment["ratio"]
        segments.append((chr, start, end, ratio))
    return segments


def get_z_score(
    segments,
    ratios,
    null_ratios,
    weights,
) -> list[float]:
    """
    Calculates between sample z-score.
    """
    zscores: list[float] = []
    for segment in segments:
        chr, start, stop, ratio = segment
        segment_null_ratio = null_ratios[chr][start:stop]
        segment_ratio = ratios[chr][start:stop]
        segment_null_ratio = [
            segment_null_ratio[i]
            for i in range(len(segment_null_ratio))
            if segment_ratio[i] != 0
        ]
        for i in range(len(segment_null_ratio)):
            for ii in range(len(segment_null_ratio[i])):
                if not np.isfinite(segment_null_ratio[i][ii]):
                    segment_null_ratio[i][ii] = np.nan
        segment_weight = weights[segment[0]][segment[1] : segment[2]]
        segment_weight = [
            segment_weight[i]
            for i in range(len(segment_weight))
            if segment_ratio[i] != 0
        ]
        null_segments = [
            np.ma.average(
                np.ma.masked_array(x, pd.isnull(x)), weights=segment_weight
            )
            for x in np.transpose(segment_null_ratio)
        ]
        null_mean = np.ma.mean([x for x in null_segments if np.isfinite(x)])
        null_sd = np.ma.std([x for x in null_segments if np.isfinite(x)])
        zscore = (segment[3] - null_mean) / null_sd
        zscore = min(zscore, 1000)
        zscore = max(zscore, -1000)
        if math.isnan(null_mean) or math.isnan(null_sd):
            zscore = "nan"
        zscores.append(zscore)
    return zscores
