# WisecondorX

import os
import re
from typing import Dict, Any, List, Tuple

import numpy as np

from wisecondorx.overall_tools import (
    exec_R,
    get_z_score,
    get_median_segment_variance,
    get_cpa,
)

"""
Writes plots.
"""


def exec_write_plots(
    outid: str,
    wd: str,
    ref_gender: str,
    beta: float,
    zscore: float,
    binsize: int,
    n_reads: int,
    cairo_flag: bool,
    ylim: str,
    regions: str,
    add_plot_title: bool,
    results: Dict[str, Any],
) -> None:
    json_plot_dir = os.path.abspath(outid + "_plot_tmp")
    json_dict = {
        "R_script": str("{}/include/plotter.R".format(wd)),
        "ref_gender": str(ref_gender),
        "beta": str(beta),
        "zscore": str(zscore),
        "binsize": str(binsize),
        "n_reads": str(n_reads),
        "cairo": str(cairo_flag),
        "results_r": results["results_r"],
        "results_w": results["results_w"],
        "results_c": results["results_c"],
        "ylim": str(ylim),
        "regions": str(regions),
        "infile": str("{}.json".format(json_plot_dir)),
        "out_dir": str("{}.plots".format(outid)),
    }

    if add_plot_title:
        # Strip away paths from the outid if need be
        json_dict["plot_title"] = str(os.path.basename(outid))

    exec_R(json_dict)


"""
Calculates zz-scores, marks aberrations and
writes tables.
"""


def generate_output_tables(
    outid: str,
    binsize: int,
    regions: str,
    bins_per_chr: list,
    ref_gender: str,
    beta: float,
    zscore: float,
    gender: str,
    n_reads: int,
    results: Dict[str, Any],
) -> None:
    _generate_bins_bed(outid, binsize, results)
    _generate_segments_and_aberrations_bed(
        outid, binsize, ref_gender, beta, zscore, results
    )
    _generate_chr_statistics_file(
        outid, bins_per_chr, binsize, gender, n_reads, results
    )
    if regions is not None:
        _generate_regions_bed(outid, binsize, regions, bins_per_chr, results)


def _generate_bins_bed(
    outid: str, binsize: int, results: Dict[str, Any]
) -> None:
    bins_file = open("{}_bins.bed".format(outid), "w")
    bins_file.write("chr\tstart\tend\tid\tratio\tzscore\n")
    results_r = results["results_r"]
    results_z = results["results_z"]
    results_z = results["results_z"]

    for chr in range(len(results_r)):
        chr_name = str(chr + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"
        feat = 1
        for i in range(len(results_r[chr])):
            r = results_r[chr][i]
            z = results_z[chr][i]
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


def _generate_regions_bed(
    outid: str,
    binsize: int,
    regions_path: str,
    bins_per_chr: List[int],
    results: Dict[str, Any],
) -> None:
    regions_file = open("{}_regions.bed".format(outid), "w")
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
            region_ratios = results["results_r"][chr][start_bin : end_bin + 1]
            region_weights = results["results_w"][chr][start_bin : end_bin + 1]
            region_zscores = results["results_z"][chr][start_bin : end_bin + 1]

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


def _generate_segments_and_aberrations_bed(
    outid: str,
    binsize: int,
    ref_gender: str,
    beta: float,
    zscore: float,
    results: Dict[str, Any],
) -> None:
    segments_file = open("{}_segments.bed".format(outid), "w")
    aberrations_file = open("{}_aberrations.bed".format(outid), "w")
    segments_file.write("chr\tstart\tend\tratio\tzscore\n")
    aberrations_file.write("chr\tstart\tend\tratio\tzscore\ttype\n")

    for segment in results["results_c"]:
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
        if (chr_name == "X" or chr_name == "Y") and ref_gender == "M":
            ploidy = 1
        if beta is not None:
            if float(segment[4]) > __get_aberration_cutoff(beta, ploidy)[1]:
                aberrations_file.write(
                    "{}\tgain\n".format("\t".join([str(x) for x in row]))
                )
            elif float(segment[4]) < __get_aberration_cutoff(beta, ploidy)[0]:
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


def __get_aberration_cutoff(beta: float, ploidy: int) -> Tuple[float, float]:
    loss_cutoff = float(np.log2((ploidy - (beta / 2)) / ploidy))
    gain_cutoff = float(np.log2((ploidy + (beta / 2)) / ploidy))
    return loss_cutoff, gain_cutoff


def _generate_chr_statistics_file(
    outid: str,
    bins_per_chr: List[int],
    binsize: int,
    gender: str,
    n_reads: int,
    results: Dict[str, Any],
) -> None:
    stats_file = open("{}_statistics.txt".format(outid), "w")
    stats_file.write("chr\tratio.mean\tratio.median\tzscore\n")
    chr_ratio_means = [
        np.ma.average(
            results["results_r"][chr], weights=results["results_w"][chr]
        )
        for chr in range(len(results["results_r"]))
    ]
    chr_ratio_medians = [
        np.median([x for x in results["results_r"][chr] if x != 0])
        for chr in range(len(results["results_r"]))
    ]

    results_c_chr = [
        [x, 0, bins_per_chr[x] - 1, chr_ratio_means[x]]
        for x in range(len(results["results_r"]))
    ]

    msv = round(
        get_median_segment_variance(
            results["results_c"], results["results_r"]
        ),
        5,
    )
    cpa = round(get_cpa(results["results_c"], binsize), 5)
    chr_z_scores = get_z_score(results_c_chr, results)

    for chr in range(len(results["results_r"])):
        chr_name = str(chr + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"

        row = [
            chr_name,
            chr_ratio_means[chr],
            chr_ratio_medians[chr],
            chr_z_scores[chr],
        ]

        stats_file.write("\t".join([str(x) for x in row]) + "\n")

    stats_file.write(
        "Gender based on --yfrac (or manually overridden by --gender): {}\n".format(
            str(gender)
        )
    )

    stats_file.write("Number of reads: {}\n".format(str(n_reads)))

    stats_file.write(
        "Standard deviation of the ratios per chromosome: {}\n".format(
            str(round(float(np.nanstd(chr_ratio_means)), 5))
        )
    )

    stats_file.write(
        "Median segment variance per bin (doi: 10.1093/nar/gky1263): {}\n".format(
            str(msv)
        )
    )

    stats_file.write(
        "Copy number profile abnormality (CPA) score (doi: 10.1186/s13073-020-00735-4): {}\n".format(
            str(cpa)
        )
    )

    stats_file.close()
