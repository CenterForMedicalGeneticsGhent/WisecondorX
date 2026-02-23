# WisecondorX

import os
import sys
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional, Union

import numpy as np
from scipy.stats import norm
import typer

import re
from wisecondorx.utils import (
    exec_R,
    get_z_score,
    scale_sample,
    gender_correct,
    get_median_segment_variance,
    get_cpa,
)
import logging


def wcx_gender(
    infile: Path = typer.Argument(..., help=".npz input file"),
    reference: Path = typer.Argument(
        ..., help="Reference .npz, as previously created with newref"
    ),
) -> None:
    """
    Returns the gender of a .npz resulting from convert, based on a Gaussian mixture model trained during the newref phase.
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
    alpha: float = typer.Option(
        1e-4, "--alpha", help="p-value cut-off for calling a CBS breakpoint."
    ),
    zscore: float = typer.Option(
        5.0, "--zscore", help="z-score cut-off for aberration calling."
    ),
    beta: Optional[float] = typer.Option(
        None,
        "--beta",
        help="When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations.",
    ),
    blacklist: Optional[str] = typer.Option(
        None,
        "--blacklist",
        help="Blacklist that masks regions in output, structure of header-less file: chr...(/t)startpos(/t)endpos(/n)",
    ),
    gender_override: Optional[str] = typer.Option(
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
    plot: bool = typer.Option(False, "--plot", help="Outputs .png plots"),
    cairo: bool = typer.Option(
        False,
        "--cairo",
        help="Uses cairo bitmap type for plotting. Might be necessary for certain setups.",
    ),
    add_plot_title: bool = typer.Option(
        False,
        "--add-plot-title",
        help="Add the output name as plot title",
    ),
    seed: Optional[int] = typer.Option(
        None, "--seed", help="Seed for segmentation algorithm"
    ),
    regions: Optional[str] = typer.Option(
        None,
        "--regions",
        help="List of regions to be marked on the output plot",
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

    if zscore <= 0:
        logging.critical(
            "Parameter --zscore should be a strictly positive number"
        )
        sys.exit(1)

    if beta is not None:
        if beta <= 0 or beta > 1:
            logging.critical(
                "Parameter --beta should be a strictly positive number lower than or equal to 1"
            )
            sys.exit(1)

    if alpha <= 0 or alpha > 1:
        logging.critical(
            "Parameter --alpha should be a strictly positive number lower than or equal to 1"
        )
        sys.exit(1)

    logging.info("Importing data ...")
    ref_file = np.load(reference, encoding="latin1", allow_pickle=True)
    sample_file = np.load(infile, encoding="latin1", allow_pickle=True)

    sample = sample_file["sample"].item()
    n_reads = sum([sum(sample[x]) for x in sample.keys()])
    ref_binsize = int(ref_file["binsize"])
    sample = scale_sample(
        sample, int(sample_file["binsize"].item()), ref_binsize
    )

    gender = predict_gender(sample, ref_file["trained_cutoff"])
    if not ref_file["is_nipt"]:
        if gender_override:
            gender = gender_override
        sample = gender_correct(sample, gender)
        ref_gender = gender
    else:
        if gender_override:
            gender = gender_override
        ref_gender = "F"

    logging.info("Normalizing autosomes ...")

    results_r, results_z, results_w, ref_sizes, m_lr, m_z = normalize(
        maskrepeats=maskrepeats,
        sample=sample,
        ref_file=ref_file,
        ref_gender="A",
    )

    if not ref_file["is_nipt"]:
        if not ref_file["has_male"] and gender == "M":
            logging.warning(
                "This sample is male, whilst the reference is created with fewer than 5 males. The female gonosomal reference will be used for X predictions."
            )
            ref_gender = "F"
        elif not ref_file["has_female"] and gender == "F":
            logging.warning(
                "This sample is female, whilst the reference is created with fewer than 5 females. The male gonosomal reference will be used for XY predictions."
            )
            ref_gender = "M"

    logging.info("Normalizing gonosomes ...")

    null_ratios_aut_per_bin = ref_file["null_ratios"]
    null_ratios_gon_per_bin = ref_file[f"null_ratios.{ref_gender}"][
        len(null_ratios_aut_per_bin) :
    ]

    results_r_2, results_z_2, results_w_2, ref_sizes_2, _, _ = normalize(
        maskrepeats=maskrepeats,
        sample=sample,
        ref_file=ref_file,
        ref_gender=ref_gender,
    )

    wd = str(os.path.dirname(os.path.realpath(__file__)))
    bpc = ref_file[f"bins_per_chr.{ref_gender}"]
    mbpc = ref_file[f"masked_bins_per_chr.{ref_gender}"]
    mbpcc = ref_file[f"masked_bins_per_chr_cum.{ref_gender}"]
    ref_mask = ref_file[f"mask.{ref_gender}"]

    rem_input: Dict[str, Any] = {
        "wd": wd,
        "binsize": ref_binsize,
        "n_reads": n_reads,
        "ref_gender": ref_gender,
        "gender": gender,
        "mask": ref_mask,
        "bins_per_chr": bpc,
        "masked_bins_per_chr": mbpc,
        "masked_bins_per_chr_cum": mbpcc,
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

    results: Dict[str, Any] = {
        "results_r": results_r,
        "results_z": results_z,
        "results_w": results_w,
        "results_nr": null_ratios,
    }

    for result_key in results.keys():
        results[result_key] = get_post_processed_result(
            minrefbins=minrefbins,
            result=results[result_key],
            ref_sizes=ref_sizes,
            rem_input=rem_input,
        )

    log_trans(results, m_lr)

    if blacklist:
        logging.info("Applying blacklist ...")
        apply_blacklist(
            blacklist_path=blacklist, binsize=ref_binsize, results=results
        )

    logging.info("Executing circular binary segmentation ...")

    results["results_c"] = exec_cbs(
        outid=outid,
        alpha=alpha,
        seed=seed,
        wd=wd,
        ref_gender=ref_gender,
        binsize=ref_binsize,
        results=results,
    )

    if bed:
        logging.info("Writing tables ...")
        generate_output_tables(
            outid=outid,
            binsize=ref_binsize,
            regions=regions,
            bins_per_chr=bpc,
            ref_gender=ref_gender,
            beta=beta,
            zscore=zscore,
            gender=gender,
            n_reads=n_reads,
            results=results,
        )

    if plot:
        logging.info("Writing plots ...")
        exec_write_plots(
            outid=outid,
            wd=wd,
            ref_gender=ref_gender,
            beta=beta,
            zscore=zscore,
            binsize=ref_binsize,
            n_reads=n_reads,
            cairo_flag=cairo,
            ylim=ylim,
            regions=regions,
            add_plot_title=add_plot_title,
            results=results,
        )

    logging.info("Finished prediction")


"""
Returns gender based on Gaussian mixture
model trained during newref phase.
"""


def predict_gender(
    sample: Dict[str, np.ndarray], trained_cutoff: float
) -> str:
    Y_fraction = float(np.sum(sample["24"])) / float(
        np.sum([np.sum(sample[x]) for x in sample.keys()])
    )
    if Y_fraction > trained_cutoff:
        return "M"
    else:
        return "F"


"""
Normalize sample for read depth and apply mask.
"""


def coverage_normalize_and_mask(
    sample: Dict[str, np.ndarray], ref_file: Dict[str, Any], ap: str
) -> np.ndarray:
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


def project_pc(
    sample_data: np.ndarray, ref_file: Dict[str, Any], ap: str
) -> np.ndarray:
    components = ref_file["pca_components{}".format(ap)]
    mean = ref_file["pca_mean{}".format(ap)]

    # Center the data
    centered = sample_data - mean

    # Project to PCA space (equivalent to pca.transform)
    transform = np.dot(centered, components.T)

    # Reconstruct from PCA space (equivalent to pca.inverse_transform)
    reconstructed = np.dot(transform, components) + mean

    return sample_data / reconstructed


"""
Defines cutoff that will add bins to a blacklist
depending on the within reference distances.
"""


def get_optimal_cutoff(ref_file: Dict[str, Any], repeats: int) -> float:
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


def normalize_repeat(
    test_data: np.ndarray,
    ref_file: Dict[str, Any],
    optimal_cutoff: float,
    ct: int,
    cp: int,
    ap: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    results_z: Optional[np.ndarray] = None
    results_r: Optional[np.ndarray] = None
    ref_sizes: Optional[np.ndarray] = None
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
    test_data: np.ndarray,
    test_copy: np.ndarray,
    ref_file: Dict[str, Any],
    optimal_cutoff: float,
    ct: int,
    cp: int,
    ap: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def get_weights(ref_file: Dict[str, Any], ap: str) -> np.ndarray:
    inverse_weights = [
        np.mean(np.sqrt(x)) for x in ref_file["distances{}".format(ap)]
    ]
    weights = np.array([1 / x for x in inverse_weights])
    return weights


"""
Unmasks results array.
"""


def inflate_results(
    results: np.ndarray, rem_input: Dict[str, Any]
) -> List[Union[float, int]]:
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


def log_trans(results: Dict[str, Any], log_r_median: float) -> None:
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


def apply_blacklist(
    blacklist_path: Optional[str], binsize: int, results: Dict[str, Any]
) -> None:
    if not blacklist_path:
        return
    blacklist = _import_bed(blacklist_path, binsize)

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


def _import_bed(
    blacklist_path: str, binsize: int
) -> Dict[int, List[List[int]]]:
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


"""
Executed CBS on results_r using results_w as weights.
Calculates segmental zz-scores.
"""


def exec_cbs(
    outid: str,
    alpha: float,
    seed: int,
    wd: str,
    ref_gender: str,
    binsize: int,
    results: Dict[str, Any],
) -> List[List[Any]]:
    json_cbs_dir = os.path.abspath(outid + "_CBS_tmp")

    json_dict = {
        "R_script": str("{}/include/CBS.R".format(wd)),
        "ref_gender": str(ref_gender),
        "alpha": str(alpha),
        "binsize": str(binsize),
        "seed": str(seed),
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


def _get_processed_cbs(
    cbs_data: List[Dict[str, Union[str, float]]],
) -> List[List[Any]]:
    results_c = []
    for i, segment in enumerate(cbs_data):
        chr = int(segment["chr"]) - 1
        s = int(segment["s"])
        e = int(segment["e"])
        r = segment["r"]
        results_c.append([chr, s, e, r])

    return results_c


"""
Control function that executes following
normalization strategies:
- coverage normalization
- between-sample normalization
- within-sample normalization
"""


def normalize(
    maskrepeats: int,
    sample: Dict[str, np.ndarray],
    ref_file: Dict[str, Any],
    ref_gender: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
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
    optimal_cutoff = get_optimal_cutoff(ref_file, maskrepeats)
    results_z, results_r, ref_sizes, m_lr, m_z = normalize_repeat(
        sample, ref_file, optimal_cutoff, ct, cp, ap
    )

    return results_r, results_z, results_w, ref_sizes, m_lr, m_z


"""
Function processes a result (e.g. results_r)
to an easy-to-interpret format. Bins without
information are set to 0.
"""


def get_post_processed_result(
    minrefbins: int,
    result: np.ndarray,
    ref_sizes: np.ndarray,
    rem_input: Dict[str, Any],
) -> List[List[Union[float, int]]]:
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
    from wisecondorx.plotter import create_plots

    plot_title = str(os.path.basename(outid)) if add_plot_title else ""
    out_dir = f"{outid}.plots"

    create_plots(
        out_dir=out_dir,
        ref_gender=ref_gender,
        beta=beta,
        zscore=zscore,
        binsize=binsize,
        n_reads=n_reads,
        cairo_flag=cairo_flag,
        ylim_str=ylim,
        regions_file=regions,
        plot_title=plot_title,
        results_r=results["results_r"],
        results_w=results["results_w"],
        results_c=results["results_c"],
    )


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
