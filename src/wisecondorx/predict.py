import typer
from pathlib import Path
import numpy as np
import logging
from wisecondorx.plotter import write_plots
import re
from scipy.stats import norm
from wisecondorx.utils import (
    sex_correct,
    scale_sample,
    Sex,
    get_cpa,
    get_z_score,
    exec_R,
    get_median_segment_variance,
)
import os
import sys
import argparse
from typing import Annotated


def wcx_sex(
    infile: Annotated[
        Path,
        typer.Argument(
            help=".npz input file",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    reference: Annotated[
        Path,
        typer.Argument(
            help="Reference .npz, as previously created with newref",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
) -> None:
    """
    Returns the sex of a .npz resulting from convert, based on a Gaussian mixture model trained during the newref phase.
    """

    ref_file = np.load(reference, encoding="latin1", allow_pickle=True)
    sample_file = np.load(infile, encoding="latin1", allow_pickle=True)
    gender = predict_sex(
        sample_file["sample"].item(), ref_file["trained_cutoff"]
    )
    if gender == "M":
        print("male")
    else:
        print("female")


def wcx_predict(
    infile: Annotated[
        Path,
        typer.Argument(
            help=".npz input file",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    reference: Annotated[
        Path,
        typer.Argument(
            help="Reference .npz, as previously created with newref",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    outid: Annotated[
        str,
        typer.Argument(
            help="Basename (w/o extension) of output files (paths are allowed, e.g. path/to/ID_1)"
        ),
    ],
    minrefbins: Annotated[
        int,
        typer.Option(
            help="Minimum amount of sensible reference bins per target bin."
        ),
    ] = 150,
    maskrepeats: Annotated[
        int,
        typer.Option(
            help="Regions with distances > mean + sd * 3 will be masked. Number of masking cycles."
        ),
    ] = 5,
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
    blacklist: Annotated[
        Path,
        typer.Option(
            help="Blacklist that masks regions in output, structure of header-less file: chr...(/t)startpos(/t)endpos(/n)"
        ),
    ] = None,
    gender: Annotated[
        Sex,
        typer.Option(
            help="Force WisecondorX to analyze this case as a male (M) or a female (F)"
        ),
    ] = None,
    ylim: Annotated[
        str, typer.Option(help="y-axis limits for plotting. e.g. [-2,2]")
    ] = "def",
    bed: Annotated[
        bool,
        typer.Option(
            help="Outputs tab-delimited .bed files, containing the most important information"
        ),
    ] = True,
    plot: Annotated[bool, typer.Option(help="Output .png plots")] = False,
    cairo: Annotated[bool, typer.Option(help="Use cairo for plotting")] = True,
    add_plot_title: Annotated[
        bool, typer.Option(help="Add the output name as plot title")
    ] = False,
    seed: Annotated[
        int, typer.Option(help="Seed for segmentation algorithm")
    ] = 42,
    regions: Annotated[
        Path,
        typer.Option(
            help="bed file with regions to be marked on the output plot"
        ),
    ] = None,
) -> None:
    """
    Find copy number aberrations.
    """
    logging.info("Starting CNA prediction")

    args = argparse.Namespace()
    args.infile = infile
    args.reference = reference
    args.outid = outid
    args.minrefbins = minrefbins
    args.maskrepeats = maskrepeats
    args.alpha = alpha
    args.zscore = zscore
    args.beta = beta
    args.blacklist = blacklist
    args.gender = gender
    args.ylim = ylim
    args.bed = bed
    args.plot = plot
    args.cairo = cairo
    args.add_plot_title = add_plot_title
    args.seed = seed
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

    if args.alpha <= 0 or args.alpha > 1:
        logging.critical(
            "Parameter --alpha should be a strictly positive number lower than or equal to 1"
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

    gender = predict_sex(sample, ref_file["trained_cutoff"])
    if not ref_file["is_nipt"]:
        if args.gender:
            gender = args.gender
        sample = sex_correct(sample, gender)
        ref_gender = gender
    else:
        if args.gender:
            gender = args.gender
        ref_gender = "F"

    logging.info("Normalizing autosomes ...")

    results_r, results_z, results_w, ref_sizes, m_lr, m_z = normalize(
        args, sample, ref_file, "A"
    )

    if not ref_file["is_nipt"]:
        if not ref_file["has_male"] and gender == "M":
            logging.warning(
                "This sample is male, whilst the reference is created with fewer than 5 males. "
                "The female gonosomal reference will be used for X predictions. Note that these might "
                "not be accurate. If the latter is desired, create a new reference and include more "
                "male samples."
            )
            ref_gender = "F"

        elif not ref_file["has_female"] and gender == "F":
            logging.warning(
                "This sample is female, whilst the reference is created with fewer than 5 females. "
                "The male gonosomal reference will be used for XY predictions. Note that these might "
                "not be accurate. If the latter is desired, create a new reference and include more "
                "female samples."
            )
            ref_gender = "M"

    logging.info("Normalizing gonosomes ...")

    null_ratios_aut_per_bin = ref_file["null_ratios"]
    null_ratios_gon_per_bin = ref_file["null_ratios.{}".format(ref_gender)][
        len(null_ratios_aut_per_bin) :
    ]

    results_r_2, results_z_2, results_w_2, ref_sizes_2, _, _ = normalize(
        args, sample, ref_file, ref_gender
    )

    rem_input = {
        "args": args,
        "wd": str(os.path.dirname(os.path.realpath(__file__))),
        "binsize": int(ref_file["binsize"]),
        "n_reads": n_reads,
        "ref_gender": ref_gender,
        "gender": gender,
        "mask": ref_file["mask.{}".format(ref_gender)],
        "bins_per_chr": ref_file["bins_per_chr.{}".format(ref_gender)],
        "masked_bins_per_chr": ref_file[
            "masked_bins_per_chr.{}".format(ref_gender)
        ],
        "masked_bins_per_chr_cum": ref_file[
            "masked_bins_per_chr_cum.{}".format(ref_gender)
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

    logging.info("Executing circular binary segmentation ...")

    results["results_c"] = exec_cbs(rem_input, results)

    if args.bed:
        logging.info("Writing tables ...")
        generate_output_tables(rem_input, results)

    if args.plot:
        logging.info("Writing plots ...")
        data = {
            "ref_gender": str(rem_input["ref_gender"]),
            "beta": str(rem_input["args"].beta),
            "zscore": str(rem_input["args"].zscore),
            "binsize": str(rem_input["binsize"]),
            "n_reads": str(rem_input["n_reads"]),
            "cairo": rem_input["args"].cairo,
            "results_r": results["results_r"],
            "results_w": results["results_w"],
            "results_c": results["results_c"],
            "ylim": str(rem_input["args"].ylim),
            "regions": str(rem_input["args"].regions),
            "out_dir": str("{}.plots".format(rem_input["args"].outid)),
        }

        if rem_input["args"].add_plot_title:
            # Strip away paths from the outid if need be
            data["plot_title"] = str(os.path.basename(rem_input["args"].outid))

        write_plots(data)

    logging.info("Finished prediction")


def normalize(args, sample, ref_file, ref_sex):
    """
    Control function that executes following
    normalization strategies:
    - coverage normalization
    - between-sample normalization
    - within-sample normalization
    """
    if ref_sex == "A":
        ap = ""
        cp = 0
        ct = 0
    else:
        ap = ".{}".format(ref_sex)
        cp = 22
        ct = ref_file["masked_bins_per_chr_cum{}".format(ap)][cp - 1]

    sample = coverage_normalize_and_mask(sample, ref_file, ap)
    sample = project_pc(sample, ref_file, ap)
    results_w = get_weights(ref_file, ap)[ct:]
    optimal_cutoff = get_optimal_cutoff(ref_file, args.maskrepeats)
    results_z, results_r, ref_sizes, m_lr, m_z = normalize_repeat(
        sample, ref_file, optimal_cutoff, ct, cp, ap
    )

    return results_r, results_z, results_w, ref_sizes, m_lr, m_z


def get_post_processed_result(args, result, ref_sizes, rem_input):
    """
    Function processes a result (e.g. results_r)
    to an easy-to-interpret format. Bins without
    information are set to 0.
    """
    infinite_mask = ref_sizes < args.minrefbins
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


def generate_output_tables(rem_input, results):
    """
    Calculates zz-scores, marks aberrations and
    writes tables.
    """
    _generate_bins_bed(rem_input, results)
    _generate_segments_and_aberrations_bed(rem_input, results)
    _generate_chr_statistics_file(rem_input, results)
    if rem_input["args"].regions is not None:
        _generate_regions_bed(rem_input, results)


def _generate_bins_bed(rem_input, results):
    bins_file = open("{}_bins.bed".format(rem_input["args"].outid), "w")
    bins_file.write("chr\tstart\tend\tid\tratio\tzscore\n")
    results_r = results["results_r"]
    results_z = results["results_z"]
    binsize = rem_input["binsize"]

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
            row = [chr_name, feat, feat + binsize - 1, feat_str, r, z]
            bins_file.write("{}\n".format("\t".join([str(x) for x in row])))
            feat += binsize
    bins_file.close()


def _generate_regions_bed(rem_input, results):
    regions_file = open("{}_regions.bed".format(rem_input["args"].outid), "w")
    regions_file.write("chr\tstart\tend\tname\tratio\tzscore\n")

    with open(rem_input["args"].regions, "r") as regions_file_handle:
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
            elif chr_name == "chrY" or chr_name == "Y":
                chr = 22
            else:
                chr = int(re.sub("chr", "", chr_name)) - 1
            start_bin = int(start) // rem_input["binsize"]
            end_bin = int(end) // rem_input["binsize"]
            if end_bin >= rem_input["bins_per_chr"][chr]:
                end_bin = rem_input["bins_per_chr"][chr] - 1

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


def _generate_segments_and_aberrations_bed(rem_input, results):
    segments_file = open(f"{rem_input['args'].outid}_segments.bed", "w")
    aberrations_file = open(f"{rem_input['args'].outid}_aberrations.bed", "w")
    segments_file.write("chr\tstart\tend\tratio\tzscore\n")
    aberrations_file.write("chr\tstart\tend\tratio\tzscore\ttype\n")

    for segment in results["results_c"]:
        chr_name = str(segment[0] + 1)
        if chr_name == "23":
            chr_name = "X"
        if chr_name == "24":
            chr_name = "Y"

        ratio_str = (
            f"{segment[4]:.3f}"
            if isinstance(segment[4], (float, np.float64, np.float32))
            else str(segment[4])
        )
        zscore_str = (
            f"{segment[3]:.3f}"
            if isinstance(segment[3], (float, np.float64, np.float32))
            else str(segment[3])
        )

        row = [
            chr_name,
            int(segment[1] * rem_input["binsize"] + 1),
            int(segment[2] * rem_input["binsize"]),
            ratio_str,
            zscore_str,
        ]
        segments_file.write(f"{'\t'.join([str(x) for x in row])}\n")

        ploidy = 2
        if (chr_name == "X" or chr_name == "Y") and rem_input[
            "ref_gender"
        ] == "M":
            ploidy = 1
        if rem_input["args"].beta is not None:
            if (
                float(segment[4])
                > _get_aberration_cutoff(rem_input["args"].beta, ploidy)[1]
            ):
                aberrations_file.write(
                    f"{'\t'.join([str(x) for x in row])}\tgain\n"
                )
            elif (
                float(segment[4])
                < _get_aberration_cutoff(rem_input["args"].beta, ploidy)[0]
            ):
                aberrations_file.write(
                    f"{'\t'.join([str(x) for x in row])}\tloss\n"
                )
        elif isinstance(segment[3], str):
            continue
        else:
            if float(segment[3]) > rem_input["args"].zscore:
                aberrations_file.write(
                    f"{'\t'.join([str(x) for x in row])}\tgain\n"
                )
            elif float(segment[3]) < -rem_input["args"].zscore:
                aberrations_file.write(
                    f"{'\t'.join([str(x) for x in row])}\tloss\n"
                )

    segments_file.close()
    aberrations_file.close()


def _get_aberration_cutoff(beta, ploidy):
    loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
    return loss_cutoff, gain_cutoff


def _generate_chr_statistics_file(rem_input, results):
    stats_file = open(f"{rem_input['args'].outid}_statistics.txt", "w")
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
            f"{chr_ratio_means[chr]:.3f}",
            f"{chr_ratio_medians[chr]:.3f}",
            f"{chr_z_scores[chr]:.3f}"
            if isinstance(chr_z_scores[chr], (float, np.float64, np.float32))
            else str(chr_z_scores[chr]),
        ]

        stats_file.write("\t".join([str(x) for x in row]) + "\n")

    stats_file.write(
        f"Gender based on --yfrac (or manually overridden by --gender): {rem_input['gender']}\n"
    )

    stats_file.write(f"Number of reads: {rem_input['n_reads']}\n")

    stats_file.write(
        f"Standard deviation of the ratios per chromosome: {np.nanstd(chr_ratio_means):.3f}\n"
    )

    stats_file.write(
        f"Median segment variance per bin (doi: 10.1093/nar/gky1263): {msv:.3f}\n"
    )

    stats_file.write(
        f"Copy number profile abnormality (CPA) score (doi: 10.1186/s13073-020-00735-4): {cpa:.3f}\n"
    )

    stats_file.close()


def predict_sex(sample, trained_cutoff):
    """
    Returns gender based on Gaussian mixture
    model trained during newref phase.
    """
    Y_fraction = float(np.sum(sample["24"])) / float(
        np.sum([np.sum(sample[x]) for x in sample.keys()])
    )
    if Y_fraction > trained_cutoff:
        return "M"
    else:
        return "F"


def coverage_normalize_and_mask(sample, ref_file, ap):
    """
    Normalize sample for read depth and apply mask.
    """
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


def project_pc(sample_data, ref_file, ap):
    """
    Project test sample to PCA space.
    """
    components = ref_file["pca_components{}".format(ap)]
    mean = ref_file["pca_mean{}".format(ap)]

    centered_data = sample_data - mean
    transform = np.dot(centered_data, components.T)

    reconstructed = np.dot(transform, components) + mean
    return sample_data / reconstructed


def get_optimal_cutoff(ref_file, repeats):
    """
    Defines cutoff that will add bins to a blacklist
    depending on the within reference distances.
    """
    distances = ref_file["distances"]
    cutoff = float("inf")
    for i in range(0, repeats):
        mask = distances < cutoff
        average = np.average(distances[mask])
        stddev = np.std(distances[mask])
        cutoff = average + 3 * stddev
    return cutoff


def normalize_repeat(test_data, ref_file, optimal_cutoff, ct, cp, ap):
    """
    Within sample normalization. Cycles through a number
    of repeats where z-scores define whether a bin is seen
    as 'normal' in a sample (in most cases this means
    'non-aberrant') -- if it is, it can be used as a
    reference for the other bins.
    """
    results_z = None
    results_r = None
    ref_sizes = None
    test_copy = np.copy(test_data)
    for i in range(3):
        results_z, results_r, ref_sizes = _normalize_once(
            test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap
        )

        test_copy[ct:][np.abs(results_z) >= norm.ppf(0.99)] = -1
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


def get_weights(ref_file, ap):
    """
    The means of sets of within-sample reference
    distances can serve as inverse weights for
    CBS, Z-scoring and plotting.
    """
    inverse_weights = np.mean(np.sqrt(ref_file[f"distances{ap}"]), axis=1)
    weights = 1.0 / inverse_weights
    # Normalize weights by their mean
    return weights / np.nanmean(weights)


def inflate_results(results, rem_input):
    """
    Unmasks results array.
    """
    mask = np.array(rem_input["mask"])
    inflated = np.zeros(len(mask))
    inflated[mask] = results
    return inflated.tolist()


def log_trans(results, log_r_median):
    """
    Log2-transforms results_r. If resulting elements are infinite,
    all corresponding possible positions (at results_r, results_z
    and results_w are set to 0 (blacklist)).
    """
    new_results_r = []
    for chr_idx in range(len(results["results_r"])):
        # Ensure it's a numpy array for vectorization
        r_orig = np.array(results["results_r"][chr_idx])
        z = np.array(results["results_z"][chr_idx])
        w = np.array(results["results_w"][chr_idx])

        # Avoid log2(0) warnings by using where or masking
        r = np.zeros_like(r_orig, dtype=float)
        valid = (r_orig > 0) & np.isfinite(r_orig)
        r[valid] = np.log2(r_orig[valid])

        # Mask for values that should be blacklisted (non-finite or non-positive original)
        to_blacklist = ~valid
        r[to_blacklist] = 0
        z[to_blacklist] = 0
        w[to_blacklist] = 0

        # Apply log_r_median to all valid (non-blacklisted) bins
        r[valid] -= log_r_median

        new_results_r.append(r.tolist())
        results["results_z"][chr_idx] = z
        results["results_w"][chr_idx] = w

    results["results_r"] = new_results_r


def apply_blacklist(rem_input, results):
    """
    Applies additional blacklist to results_r, results_z
    and results_w if requested.
    """
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


def exec_cbs(rem_input, results):
    """
    Executed CBS on results_r using results_w as weights.
    Calculates segmental zz-scores.
    """
    json_cbs_dir = os.path.abspath(rem_input["args"].outid + "_CBS_tmp")

    json_dict = {
        "R_script": str("{}/include/CBS.R".format(rem_input["wd"])),
        "ref_gender": str(rem_input["ref_gender"]),
        "alpha": str(rem_input["args"].alpha),
        "binsize": str(rem_input["binsize"]),
        "seed": str(rem_input["args"].seed),
        "results_r": results["results_r"],
        "results_w": results["results_w"],
        "infile": str("{}_01.json".format(json_cbs_dir)),
        "outfile": str("{}_02.json".format(json_cbs_dir)),
    }

    results_c = _get_processed_cbs(exec_R(json_dict))
    segment_z = get_z_score(results_c, results)
    results_c = [
        results_c[i][:3] + [segment_z[i]] + [results_c[i][3]]
        for i in range(len(results_c))
    ]
    return results_c


def _get_processed_cbs(cbs_data):
    results_c = []
    for i, segment in enumerate(cbs_data):
        chr = int(segment["chr"]) - 1
        s = int(segment["s"])
        e = int(segment["e"])
        r = segment["r"]
        results_c.append([chr, s, e, r])

    return results_c
