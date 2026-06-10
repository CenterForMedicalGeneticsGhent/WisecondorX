# WisecondorX

import re
import numpy as np

from wisecondorx.overall_tools import (
    get_z_score,
    get_median_segment_variance,
    get_cpa,
)


"""
Calculates zz-scores, marks aberrations and
writes tables.
"""


def generate_output_tables(cna_df, segments_df, args):
    _generate_bins_bed(cna_df, args)
    _generate_segments_bed(segments_df, args)
    _generate_aberrations_bed(segments_df, args)
    # _generate_chr_statistics_file(segments_df, args)
    # if args.regions is not None:
    #     _generate_regions_bed(segments_df, args)


def _generate_bins_bed(cna_df, args):
    with open(f"{args.outid}_bins.bed", "w") as bins_file:
        cna_df.to_csv(
            bins_file,
            sep="\t",
            index=False,
            columns=["chrom", "start", "end", "ratios", "zscores"],
            header=["chr", "start", "end", "ratio", "zscore"],
            mode="w",
        )


def _generate_regions_bed(cna_df, segments_df, args):
    regions_file = open("{}_regions.bed".format(args.outid), "w")
    regions_file.write("chr\tstart\tend\tname\tratio\tzscore\n")

    with open(args.regions, "r") as regions_file_handle:
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
            start_bin = int(start) // args.binsize
            end_bin = int(end) // args.binsize
            if end_bin >= args.bins_per_chr[chr]:
                end_bin = args.bins_per_chr[chr] - 1

            if start_bin < 0 or end_bin < 0 or start_bin > end_bin:
                regions_file.write(
                    "Skipping invalid region: {}\n".format("\t".join(region))
                )
                continue

            # Extract ratios, weights, and z-scores for the region
            region_ratios = cna_df[chr][start_bin : end_bin + 1]
            region_weights = cna_df[chr][start_bin : end_bin + 1]
            region_zscores = cna_df[chr][start_bin : end_bin + 1]

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


def _generate_segments_bed(segments_df, args):
    with open(f"{args.outid}_segments.bed", "w") as segments_file:
        segments_df.to_csv(
            segments_file,
            sep="\t",
            index=False,
            columns=["chrom", "start", "end", "ratio", "zscore"],
            header=["chr", "start", "end", "ratio", "zscore"],
            mode="w",
        )


def _generate_aberrations_bed(segments_df, args):
    aberrations_df = segments_df.copy()

    chrom_norm = aberrations_df["chrom"].astype(str)
    is_gonosomal = chrom_norm.isin(["X", "Y", "23", "24"])

    gender_value = getattr(args.gender, "value", args.gender)
    is_male = str(gender_value).upper() == "M"
    ploidy = np.where(is_male & is_gonosomal, 1.0, 2.0)

    aberrations_df["type"] = "neutral"
    if args.beta is not None:
        loss_cutoff = np.log2((ploidy - (args.beta / 2.0)) / ploidy)
        gain_cutoff = np.log2((ploidy + (args.beta / 2.0)) / ploidy)
        aberrations_df.loc[aberrations_df["ratio"] > gain_cutoff, "type"] = (
            "gain"
        )
        aberrations_df.loc[aberrations_df["ratio"] < loss_cutoff, "type"] = (
            "loss"
        )
    else:
        aberrations_df = aberrations_df[aberrations_df["zscore"].notna()]
        aberrations_df.loc[aberrations_df["zscore"] > args.zscore, "type"] = (
            "gain"
        )
        aberrations_df.loc[aberrations_df["zscore"] < -args.zscore, "type"] = (
            "loss"
        )

    aberrations_df = aberrations_df[
        aberrations_df["type"].isin(["gain", "loss"])
    ]

    with open(f"{args.outid}_aberrations.bed", "w") as aberrations_file:
        aberrations_df.to_csv(
            aberrations_file,
            sep="\t",
            index=False,
            columns=["chrom", "start", "end", "ratio", "zscore", "type"],
            header=["chr", "start", "end", "ratio", "zscore", "type"],
            mode="w",
        )


def _generate_chr_statistics_file(rem_input, results):
    stats_file = open("{}_statistics.txt".format(rem_input["args"].outid), "w")
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
        [x, 0, rem_input["bins_per_chr"][x] - 1, chr_ratio_means[x]]
        for x in range(len(results["results_r"]))
    ]

    msv = round(
        get_median_segment_variance(
            results["results_c"], results["results_r"]
        ),
        5,
    )
    cpa = round(get_cpa(results["results_c"], rem_input["binsize"]), 5)
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
            str(rem_input["gender"])
        )
    )

    stats_file.write("Number of reads: {}\n".format(str(rem_input["n_reads"])))

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
