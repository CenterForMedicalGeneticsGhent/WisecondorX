# WisecondorX

import re
import numpy as np
import os

from scipy.stats import norm
from pathlib import Path
from tempfile import TemporaryDirectory
import logging
import json
import subprocess
import sys
import typer
from typing import Annotated
import argparse

from wisecondorx.overall_tools import (
    gender_correct,
    get_z_score,
    get_median_segment_variance,
    scale_sample,
)
from wisecondorx.plotter import write_plots
from wisecondorx.overall_tools import Sex


def wcx_sex(
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
    gender = predict_sex(
        sample_file["sample"].item(), ref_file["trained_cutoff"]
    )
    if gender == "M":
        print("male")
    else:
        print("female")


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
    alpha: Annotated[
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
            "--cbs-nperm",
            help="Number of permutations used by CBS.",
            min=1,
        ),
    ] = 10000,
    cbs_min_width: Annotated[
        int,
        typer.Option(
            "--cbs-min-width",
            help="Minimum markers for a CBS changed segment (2-5).",
            min=2,
            max=5,
        ),
    ] = 2,
    cbs_undo_splits: str = typer.Option(
        "none",
        "--cbs-undo-splits",
        help="Undo CBS splits using one of: none, prune, sdundo.",
    ),
    cbs_undo_sd: Annotated[
        float,
        typer.Option(
            "--cbs-undo-sd",
            help="SD cutoff used when --cbs-undo-splits=sdundo.",
            min=0.0,
        ),
    ] = 3.0,
    cbs_smooth: bool = typer.Option(
        False,
        "--cbs-smooth",
        help="Apply DNAcopy smooth.CNA before segmentation.",
    ),
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

    args = argparse.Namespace()
    args.infile = infile
    args.reference = reference
    args.outid = outid
    args.minrefbins = minrefbins
    args.maskrepeats = maskrepeats
    args.alpha = alpha
    args.cbs_nperm = cbs_nperm
    args.cbs_min_width = cbs_min_width
    args.cbs_undo_splits = cbs_undo_splits
    args.cbs_undo_sd = cbs_undo_sd
    args.cbs_smooth = cbs_smooth
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

    if args.cbs_undo_splits not in {"none", "prune", "sdundo"}:
        logging.critical(
            "Parameter --cbs-undo-splits should be one of: none, prune, sdundo"
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

    gender = Sex(predict_sex(sample, ref_file["trained_cutoff"]))
    if not ref_file["is_nipt"]:
        if args.gender:
            gender = args.gender
        sample = gender_correct(sample, gender)
        ref_gender = gender.value
    else:
        if args.gender:
            gender = args.gender
        ref_gender = "F"

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
            ref_gender = "F"

        elif not ref_file["has_female"] and gender == Sex.FEMALE:
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


def predict_sex(sample, trained_cutoff):
    """
    Returns sex based on Gaussian mixture
    model trained during newref phase.
    """
    Y_fraction = float(np.sum(sample["24"])) / float(
        np.sum([np.sum(sample[x]) for x in sample.keys()])
    )
    return "M" if Y_fraction > trained_cutoff else "F"


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
    for _ in range(repeats):
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
    test_copy = np.copy(test_data)
    results_z, results_r, ref_sizes = _normalize_once(
        test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap
    )
    test_copy[ct:][np.abs(results_z) >= norm.ppf(0.99)] = -1

    for _ in range(2):
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
    inverse_weights = [
        np.mean(np.sqrt(x)) for x in ref_file["distances{}".format(ap)]
    ]
    weights = np.array([1 / x for x in inverse_weights])
    return weights


def inflate_results(results, rem_input):
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


def log_trans(results, log_r_median):
    """
    Log2-transforms results_r. If resulting elements are infinite,
    all corresponding possible positions (at results_r, results_z
    and results_w are set to 0 (blacklist)).
    """

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


def apply_blacklist(rem_input, results):
    """
    Applies additional blacklist to results_r, results_z
    and results_w if requested.
    """
    blacklist = _import_bed(rem_input)
    results_r = results["results_r"]
    results_z = results["results_z"]
    results_w = results["results_w"]
    has_y = len(results_r) >= 24

    for chrom, spans in blacklist.items():
        for start, end in spans:
            for pos in range(start, end):
                if not has_y and chrom == 23:
                    continue
                if pos >= len(results_r[chrom]) or pos < 0:
                    continue
                results_r[chrom][pos] = 0
                results_z[chrom][pos] = 0
                results_w[chrom][pos] = 0


def _import_bed(rem_input):
    bed = {}
    binsize = rem_input["binsize"]
    with open(rem_input["args"].blacklist) as bed_handle:
        for line in bed_handle:
            chr_name, s, e = line.strip().split("\t")
            if chr_name.startswith("chr"):
                chr_name = chr_name[3:]
            if chr_name == "X":
                chr_name = "23"
            if chr_name == "Y":
                chr_name = "24"

            chrom = int(chr_name) - 1
            bed.setdefault(chrom, []).append(
                [
                    int(int(s) / binsize),
                    int(int(e) / binsize) + 1,
                ]
            )
    return bed


def exec_cbs(rem_input, results):
    """
    Executed CBS on results_r using results_w as weights.
    Calculates segmental zz-scores.
    """
    script_path = Path(__file__).parent / "include" / "CBS.R"
    sample_id = getattr(rem_input["args"], "outid", "sample")
    results_c = []
    with TemporaryDirectory() as tmpdir:
        cbs_input = Path(tmpdir) / "cbs_input.json"
        cbs_output = Path(tmpdir) / "cbs_output.json"

        json.dump(
            {
                "ratios": results["results_r"],
                "weights": results["results_w"],
            },
            cbs_input.open("w"),
            indent=4,
        )
        r_cmd = [
            "Rscript",
            str(script_path),
            "--infile",
            str(cbs_input),
            "--outfile",
            str(cbs_output),
            "--id",
            str(os.path.basename(sample_id)),
            "--alpha",
            str(rem_input["args"].alpha),
            "--cbs_nperm",
            str(rem_input["args"].cbs_nperm),
            "--cbs_min_width",
            str(rem_input["args"].cbs_min_width),
            "--cbs_undo_splits",
            str(rem_input["args"].cbs_undo_splits),
            "--cbs_undo_sd",
            str(rem_input["args"].cbs_undo_sd),
            "--cbs_smooth",
            str(rem_input["args"].cbs_smooth).upper(),
            "--seed",
            str(rem_input["args"].seed),
        ]
        logging.debug("CBS cmd: {}".format(r_cmd))

        try:
            subprocess.check_call(r_cmd)
            with cbs_output.open("r") as f:
                results_c = _get_processed_cbs(json.load(f))
        except subprocess.CalledProcessError as e:
            logging.critical(f"Rscript failed: {e}")
            cbs_input.rename(os.getcwd() + "/cbs_input.json")
            sys.exit(1)

    segment_z = get_z_score(results_c, results)
    results_c = [
        segment[:3] + [z] + [segment[3]]
        for segment, z in zip(results_c, segment_z)
    ]
    return results_c  # list of [chr, start_bin, end_bin, segment_z_score, segment_mean_ratio]


def _get_processed_cbs(cbs_data):
    results_c = []
    for segment in cbs_data:
        chrom = int(segment["chrom"]) - 1
        start = int(segment["loc.start"])
        end = int(segment["loc.end"])
        seg_mean = float(segment["seg.mean"])
        results_c.append([chrom, start, end, seg_mean])

    return results_c


def normalize(args, sample, ref_file, ref_gender):
    """
    Control function that executes following
    normalization strategies:
    - coverage normalization
    - between-sample normalization
    - within-sample normalization
    """
    if ref_gender == "A":
        ap = ""
        cp = 0
        ct = 0
    else:
        ap = ".{}".format(ref_gender)
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
    with open("{}_bins.bed".format(rem_input["args"].outid), "w") as bins_file:
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


def _generate_regions_bed(rem_input, results):
    with open("{}_regions.bed".format(rem_input["args"].outid), "w") as regions_file:
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
            if chr_name == "chrY" or chr_name == "Y":
                chr = 22
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



def _generate_segments_and_aberrations_bed(rem_input, results):
    with open(
        "{}_segments.bed".format(rem_input["args"].outid), "w"
    ) as segments_file, open(
        "{}_aberrations.bed".format(rem_input["args"].outid), "w"
    ) as aberrations_file:
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
                int(segment[1] * rem_input["binsize"] + 1),
                int(segment[2] * rem_input["binsize"]),
                segment[4],
                segment[3],
            ]
            segments_file.write("{}\n".format("\t".join([str(x) for x in row])))

            ploidy = 2
            if (chr_name == "X" or chr_name == "Y") and rem_input[
                "ref_gender"
            ] == "M":
                ploidy = 1
            if rem_input["args"].beta is not None:
                if (
                    float(segment[4])
                    > __get_aberration_cutoff(rem_input["args"].beta, ploidy)[1]
                ):
                    aberrations_file.write(
                        "{}\tgain\n".format("\t".join([str(x) for x in row]))
                    )
                elif (
                    float(segment[4])
                    < __get_aberration_cutoff(rem_input["args"].beta, ploidy)[0]
                ):
                    aberrations_file.write(
                        "{}\tloss\n".format("\t".join([str(x) for x in row]))
                    )
            elif isinstance(segment[3], str):
                continue
            else:
                if float(segment[3]) > rem_input["args"].zscore:
                    aberrations_file.write(
                        "{}\tgain\n".format("\t".join([str(x) for x in row]))
                    )
                elif float(segment[3]) < -rem_input["args"].zscore:
                    aberrations_file.write(
                        "{}\tloss\n".format("\t".join([str(x) for x in row]))
                    )


def __get_aberration_cutoff(beta, ploidy):
    loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
    return loss_cutoff, gain_cutoff


def _generate_chr_statistics_file(rem_input, results):
    with open("{}_statistics.txt".format(rem_input["args"].outid), "w") as stats_file:
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
        cpa = round(_get_cpa(results["results_c"], rem_input["binsize"]), 5)
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


def _get_cpa(results_c, binsize):
    """
    Returns CPA, measure for sample-wise abnormality.
    """
    x = 0
    for _, start, stop, _, seg_mean in results_c:
        x += (stop - start + 1) * binsize * abs(seg_mean)
    CPA = x / len(results_c) * (10**-8)
    return CPA
