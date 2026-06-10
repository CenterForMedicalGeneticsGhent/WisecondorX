#!/usr/bin/env python3

import logging
import os
import sys
import warnings
from typing import Annotated
import numpy as np
import pandas as pd
import typer
import argparse
from pathlib import Path
from wisecondorx.convert_tools import wcx_convert
from wisecondorx.newref_control import (
    tool_newref_prep,
    tool_newref_main,
    tool_newref_merge,
)
from wisecondorx.newref_tools import train_gender_model, get_mask
from wisecondorx.overall_tools import gender_correct, scale_sample, Sex
from wisecondorx.predict_control import normalize, get_post_processed_result
from wisecondorx.predict_tools import (
    get_segment_zscore,
    log_trans,
    exec_cbs,
    apply_blacklist,
    predict_gender,
)
from wisecondorx.predict_output import generate_output_tables
from wisecondorx.plotter import write_plots
from wisecondorx.ref_qc import qc_reference

from wisecondorx import __version__

VERSION: str = __version__


def setup_logging(loglevel: str = "INFO") -> None:
    logging.basicConfig(
        format="[%(levelname)s - %(asctime)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=getattr(logging, loglevel.upper(), None),
    )


def _build_branch_masks(samples, genders, is_nipt):
    total_mask, bins_per_chr = get_mask(samples)

    bins_per_chr = np.array(bins_per_chr, dtype=int)
    chr22_end = int(np.sum(bins_per_chr[:22]))
    chr23_end = int(np.sum(bins_per_chr[:23]))
    chr24_end = int(np.sum(bins_per_chr[:24]))

    mask_a = np.zeros_like(total_mask, dtype=bool)
    mask_a[:chr22_end] = total_mask[:chr22_end]

    mask_f_sex = None
    if genders.count("F") > 4:
        female_mask, _ = get_mask(samples[np.array(genders) == "F"])
        mask_f_sex = np.zeros_like(total_mask, dtype=bool)
        if len(female_mask) >= chr23_end:
            mask_f_sex[chr22_end:chr23_end] = female_mask[chr22_end:chr23_end]

    mask_m_sex = None
    if genders.count("M") > 4 and not is_nipt:
        male_mask, _ = get_mask(samples[np.array(genders) == "M"])
        mask_m_sex = np.zeros_like(total_mask, dtype=bool)
        if len(male_mask) >= chr23_end:
            mask_m_sex[chr22_end:chr23_end] = male_mask[chr22_end:chr23_end]
        if len(male_mask) >= chr24_end:
            mask_m_sex[chr23_end:chr24_end] = male_mask[chr23_end:chr24_end]

    return bins_per_chr.tolist(), mask_a, mask_f_sex, mask_m_sex


def wcx_newref(
    infiles: list[Path] = typer.Argument(
        ...,
        help="Path to all reference data files (e.g. path/to/reference/*.npz)",
    ),
    outfile: Path = typer.Argument(
        ...,
        help="Path and filename for the reference output (e.g. path/to/myref.npz)",
    ),
    nipt: bool = typer.Option(False, "--nipt", help="Use flag for NIPT"),
    yfrac: Annotated[
        float,
        typer.Option(
            "--yfrac",
            min=0.0,
            max=1.0,
            help="Use to manually set the Y read fraction cutoff, which defines sex",
        ),
    ] = None,
    plotyfrac: Path = typer.Option(
        None,
        "--plotyfrac",
        help="Path to yfrac .png plot for optimization; software will stop after plotting",
    ),
    refsize: int = typer.Option(
        300, "--refsize", help="Amount of reference locations per target"
    ),
    binsize: int = typer.Option(
        5000, "--binsize", help="Size of target bins in base pairs"
    ),
    cpus: int = typer.Option(
        1, "--cpus", help="Use multiple cores to find reference bins"
    ),
) -> None:
    """
    Create a new reference using healthy reference samples.
    """
    logging.info("Creating new reference")

    args = argparse.Namespace()
    args.infiles = infiles
    args.outfile = outfile
    args.nipt = nipt
    args.yfrac = yfrac
    args.plotyfrac = plotyfrac
    args.refsize = refsize
    args.binsize = binsize
    args.cpus = cpus

    if args.yfrac is not None:
        if args.yfrac < 0 or args.yfrac > 1:
            logging.critical(
                "Parameter --yfrac should be a positive number lower than or equal to 1"
            )
            sys.exit()

    split_path = list(os.path.split(args.outfile))
    if split_path[-1][-4:] == ".npz":
        split_path[-1] = split_path[-1][:-4]
    base_path = os.path.join(split_path[0], split_path[1])

    args.basepath = base_path
    args.prepfile = "{}_prep.npz".format(base_path)
    args.prepdatafile = "{}_prep_data.npy".format(base_path)
    args.partfile = "{}_part".format(base_path)

    samples = []
    logging.info("Importing data ...")
    for infile in args.infiles:
        logging.info("Loading: {}".format(infile))
        npzdata = np.load(infile, encoding="latin1", allow_pickle=True)
        sample = npzdata["sample"].item()
        binsize = int(npzdata["binsize"])
        logging.info("Binsize: {}".format(int(binsize)))
        samples.append(scale_sample(sample, binsize, args.binsize))

    samples = np.array(samples)
    genders, trained_cutoff = train_gender_model(args, samples)

    if genders.count("F") < 5 and args.nipt:
        logging.warning(
            "A NIPT reference should have at least 5 female feti samples. Removing --nipt flag."
        )
        args.nipt = False
    if not args.nipt:
        for i, sample in enumerate(samples):
            samples[i] = gender_correct(sample, genders[i])

    bins_per_chr, mask_A, mask_F_sex, mask_M_sex = _build_branch_masks(
        samples, genders, args.nipt
    )
    chr22_end = int(np.sum(bins_per_chr[:22]))
    chr23_end = int(np.sum(bins_per_chr[:23]))
    chr24_end = int(np.sum(bins_per_chr[:24]))

    outfiles = []
    if len(genders) > 9:
        logging.info("Starting autosomal reference creation ...")
        args.tmpoutfile = "{}.tmp.A.npz".format(args.basepath)
        outfiles.append(args.tmpoutfile)
        tool_newref_prep(args, samples, "A", mask_A, bins_per_chr)
        logging.info("This might take a while ...")
        tool_newref_main(args, args.cpus)
    else:
        logging.critical(
            "Provide at least 10 samples to enable the generation of a reference."
        )
        sys.exit()

    if genders.count("F") > 4:
        logging.info("Starting female gonosomal reference creation ...")
        args.tmpoutfile = "{}.tmp.F.npz".format(args.basepath)
        outfiles.append(args.tmpoutfile)
        mask_F = np.copy(mask_A)
        if mask_F_sex is not None:
            mask_F[chr22_end:chr23_end] = mask_F_sex[chr22_end:chr23_end]
        tool_newref_prep(
            args, samples[np.array(genders) == "F"], "F", mask_F, bins_per_chr
        )
        logging.info("This might take a while ...")
        tool_newref_main(args, 1)
    else:
        logging.warning(
            "Provide at least 5 female samples to enable normalization of female gonosomes."
        )

    if not args.nipt:
        if genders.count("M") > 4:
            logging.info("Starting male gonosomal reference creation ...")
            args.tmpoutfile = "{}.tmp.M.npz".format(args.basepath)
            outfiles.append(args.tmpoutfile)
            mask_M = np.copy(mask_A)
            if mask_M_sex is not None:
                mask_M[chr22_end:chr24_end] = mask_M_sex[chr22_end:chr24_end]
            tool_newref_prep(
                args,
                samples[np.array(genders) == "M"],
                "M",
                mask_M,
                bins_per_chr,
            )
            tool_newref_main(args, 1)
        else:
            logging.warning(
                "Provide at least 5 male samples to enable normalization of male gonosomes."
            )

    tool_newref_merge(args, outfiles, trained_cutoff)

    logging.info("Running QC on the newly created reference...")
    qc_reference(args.outfile)

    logging.info("Finished creating reference")


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


app = typer.Typer(
    name="wisecondorx",
    help="WisecondorX: Copy Number Aberration detection from Whole Genome Sequencing data.",
    add_completion=False,
)
app.command(name="convert")(wcx_convert)
app.command(name="newref")(wcx_newref)
app.command(name="gender")(wcx_gender)
app.command(name="predict")(wcx_predict)


@app.callback()
def main_callback(
    loglevel: str = typer.Option(
        "INFO",
        "--loglevel",
        help="Logging level (info, warning, debug, error, critical)",
    ),
) -> None:
    warnings.filterwarnings("ignore")
    setup_logging(loglevel=loglevel)


if __name__ == "__main__":
    app()
