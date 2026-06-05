# WisecondorX

import json
import logging
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
from wisecondorx.overall_tools import get_z_score

"""
Returns gender based on Gaussian mixture
model trained during newref phase.
"""


def predict_gender(sample, trained_cutoff):
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


def exec_cbs(rem_input, results):
    """
    Executed CBS on results_r using results_w as weights.
    Calculates segmental zz-scores.
    """

    # script params
    script_path = Path(__file__).parent / "include" / "CBS.R"
    sample_id = getattr(rem_input["args"], "outid", "sample")

    # prepare dataframe containing vectors for R script
    # convert 0 weights to very small numbers to avoid issues with CBS
    # remove chromosomes with only 0 ratios
    # chrom     :  [Integer Vector]  -> e.g., 1, 1, 1, 2, 2, 3...
    # maploc    :  [Numeric Vector]  -> e.g., 1000, 2000, 3000, 1000, 2000...
    # genomedat  :  [Numeric Vector]  -> e.g., -0.05, 0.12, -0.44, 0.02...
    # weights   :  [Numeric Vector]  -> e.g., 0.5, 0.8, 0.3, 0.9, 0.7...
    chrom_names = list(range(1, len(results["results_r"]) + 1))
    cna_df = (
        pd.DataFrame(
            [
                {
                    "chrom": chrom,
                    "maploc": (idx + 1) * rem_input["binsize"],
                    "genomedat": float(ratio),
                    "weights": np.nextafter(0, 1)
                    if weight == 0
                    else float(weight),
                }
                for chrom, ratios, weights in zip(
                    chrom_names, results["results_r"], results["results_w"]
                )
                for idx, (ratio, weight) in enumerate(zip(ratios, weights))
            ]
        )
        .groupby("chrom")
        .filter(lambda x: x["genomedat"].notna().any())
        .sort_values(by=["chrom", "maploc"], ascending=[True, True])
    )

    # Run CBS and process results
    results_c = []
    with TemporaryDirectory() as tmpdir:
        cbs_input = Path(tmpdir) / "cbs_input.json"
        cbs_output = Path(tmpdir) / "cbs_output.json"

        with cbs_input.open("w") as f:
            json.dump(
                cna_df.replace({np.nan: None}).to_dict(orient="list"),
                f,
                indent=4,
            )

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
            str(os.path.basename(sample_id)),
            "--seed",
            str(rem_input["args"].cbs_seed),
            # CBS
            "--alpha",
            str(rem_input["args"].cbs_alpha),
            "--nperm",
            str(rem_input["args"].cbs_nperm),
            "--p_method",
            str(rem_input["args"].cbs_p_method),
            "--min_width",
            str(rem_input["args"].cbs_min_width),
            "--kmax",
            str(rem_input["args"].cbs_kmax),
            "--nmin",
            str(rem_input["args"].cbs_nmin),
            "--eta",
            str(rem_input["args"].cbs_eta),
            "--trim",
            str(rem_input["args"].cbs_trim),
            # set verbose level based on Python logging level (debug -> 3, else 0)
            "--verbose",
            "3" if logging.getLogger().getEffectiveLevel() == 10 else "0",
        ]
        logging.debug("CBS cmd: {}".format(r_cmd))

        try:
            subprocess.check_call(r_cmd)
            with cbs_output.open("r") as f:
                results_c = pd.read_json(f)
        except subprocess.CalledProcessError as e:
            logging.critical(f"Rscript failed: {e}")
            cbs_input.rename(os.getcwd() + "/cbs_input.json")
            sys.exit(1)

    segment_z = get_z_score(results_c, results)
    results_c = [
        results_c[i][:3] + [segment_z[i]] + [results_c[i][3]]
        for i in range(len(results_c))
    ]
    return results_c
