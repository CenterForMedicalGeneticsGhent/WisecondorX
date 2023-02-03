#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import warnings

import numpy as np

from wisecondorx.convert_tools import convert_reads
from wisecondorx.newref_control import (
    tool_newref_prep,
    tool_newref_main,
    tool_newref_merge,
)
from wisecondorx.newref_tools import train_gender_model, get_mask
from wisecondorx.overall_tools import gender_correct, scale_sample
from wisecondorx.predict_control import normalize, get_post_processed_result
from wisecondorx.predict_output import generate_output_tables, exec_write_plots
from wisecondorx.predict_tools import (
    log_trans,
    exec_cbs,
    apply_blacklist,
    predict_gender,
)


def tool_convert(args):
    logging.info("Starting conversion")

    sample, qual_info = convert_reads(args)
    np.savez_compressed(
        args.outfile, binsize=args.binsize, sample=sample, quality=qual_info
    )

    logging.info("Finished conversion")


def tool_newref(args):
    logging.info("Creating new reference")

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

    total_mask, bins_per_chr = get_mask(samples)
    if genders.count("F") > 4:
        mask_F, _ = get_mask(samples[np.array(genders) == "F"])
        total_mask = total_mask & mask_F
    if genders.count("M") > 4 and not args.nipt:
        mask_M, _ = get_mask(samples[np.array(genders) == "M"])
        total_mask = total_mask & mask_M

    outfiles = []
    if len(genders) > 9:
        logging.info("Starting autosomal reference creation ...")
        args.tmpoutfile = "{}.tmp.A.npz".format(args.basepath)
        outfiles.append(args.tmpoutfile)
        tool_newref_prep(args, samples, "A", total_mask, bins_per_chr)
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
        tool_newref_prep(
            args, samples[np.array(genders) == "F"], "F", total_mask, bins_per_chr
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
            tool_newref_prep(
                args, samples[np.array(genders) == "M"], "M", total_mask, bins_per_chr
            )
            tool_newref_main(args, 1)
        else:
            logging.warning(
                "Provide at least 5 male samples to enable normalization of male gonosomes."
            )

    tool_newref_merge(args, outfiles, trained_cutoff)

    logging.info("Finished creating reference")


def tool_test(args):
    logging.info("Starting CNA prediction")

    if not args.bed and not args.plot:
        logging.critical(
            "No output format selected. "
            "Select at least one of the supported output formats (--bed, --plot)"
        )
        sys.exit()

    if args.zscore <= 0:
        logging.critical("Parameter --zscore should be a strictly positive number")
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

    gender = predict_gender(sample, ref_file["trained_cutoff"])
    if not ref_file["is_nipt"]:
        if args.gender:
            gender = args.gender
        sample = gender_correct(sample, gender)
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
        "masked_bins_per_chr": ref_file["masked_bins_per_chr.{}".format(ref_gender)],
        "masked_bins_per_chr_cum": ref_file[
            "masked_bins_per_chr_cum.{}".format(ref_gender)
        ],
    }

    del ref_file

    results_r = np.append(results_r, results_r_2)
    results_z = np.append(results_z, results_z_2) - m_z
    results_w = np.append(
        results_w * np.nanmean(results_w_2), results_w_2 * np.nanmean(results_w)
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
        + [x.tolist() for x in null_ratios_gon_per_bin]
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
        exec_write_plots(rem_input, results)

    logging.info("Finished prediction")


def output_gender(args):
    ref_file = np.load(args.reference, encoding="latin1", allow_pickle=True)
    sample_file = np.load(args.infile, encoding="latin1", allow_pickle=True)
    gender = predict_gender(sample_file["sample"].item(), ref_file["trained_cutoff"])
    if gender == "M":
        print("male")
    else:
        print("female")


def main():
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser(description="WisecondorX")
    parser.add_argument(
        "--loglevel",
        type=str,
        default="INFO",
        choices=["info", "warning", "debug", "error", "critical"],
    )
    subparsers = parser.add_subparsers()

    parser_convert = subparsers.add_parser(
        "convert",
        description="Convert and filter a aligned reads to .npz",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_convert.add_argument(
        "infile", type=str, help="aligned reads input for conversion"
    )
    parser_convert.add_argument("outfile", type=str, help="Output .npz file")
    parser_convert.add_argument(
        "-r",
        "--reference",
        type=str,
        help="Fasta reference to be used during cram conversion",
    )
    parser_convert.add_argument(
        "--binsize", type=float, default=5e3, help="Bin size (bp)"
    )
    parser_convert.add_argument(
        "--normdup", action="store_true", help="Do not remove duplicates"
    )
    parser_convert.set_defaults(func=tool_convert)

    parser_newref = subparsers.add_parser(
        "newref",
        description="Create a new reference using healthy reference samples",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_newref.add_argument(
        "infiles",
        type=str,
        nargs="+",
        help="Path to all reference data files (e.g. path/to/reference/*.npz)",
    )
    parser_newref.add_argument(
        "outfile",
        type=str,
        help="Path and filename for the reference output (e.g. path/to/myref.npz)",
    )
    parser_newref.add_argument(
        "--nipt",
        action="store_true",
        help="Use flag for NIPT (e.g. path/to/reference/*.npz)",
    )
    parser_newref.add_argument(
        "--yfrac",
        type=float,
        default=None,
        help="Use to manually set the y read fraction cutoff, which defines gender",
    )
    parser_newref.add_argument(
        "--plotyfrac",
        type=str,
        default=None,
        help="Path to yfrac .png plot for --yfrac optimization (e.g. path/to/plot.png); software will stop after plotting after which --yfrac can be set manually",
    )
    parser_newref.add_argument(
        "--refsize",
        type=int,
        default=300,
        help="Amount of reference locations per target",
    )
    parser_newref.add_argument(
        "--binsize",
        type=int,
        default=1e5,
        help="Scale samples to this binsize, multiples of existing binsize only",
    )
    parser_newref.add_argument(
        "--cpus", type=int, default=1, help="Use multiple cores to find reference bins"
    )
    parser_newref.set_defaults(func=tool_newref)

    parser_gender = subparsers.add_parser(
        "gender",
        description="Returns the gender of a .npz resulting from convert, "
        "based on a Gaussian mixture model trained during the newref phase",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_gender.add_argument("infile", type=str, help=".npz input file")
    parser_gender.add_argument(
        "reference", type=str, help="Reference .npz, as previously created with newref"
    )
    parser_gender.set_defaults(func=output_gender)

    parser_test = subparsers.add_parser(
        "predict",
        description="Find copy number aberrations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_test.add_argument("infile", type=str, help=".npz input file")
    parser_test.add_argument(
        "reference", type=str, help="Reference .npz, as previously created with newref"
    )
    parser_test.add_argument(
        "outid",
        type=str,
        help="Basename (w/o extension) of output files (paths are allowed, e.g. path/to/ID_1)",
    )
    parser_test.add_argument(
        "--minrefbins",
        type=int,
        default=150,
        help="Minimum amount of sensible reference bins per target bin.",
    )
    parser_test.add_argument(
        "--maskrepeats",
        type=int,
        default=5,
        help="Regions with distances > mean + sd * 3 will be masked. Number of masking cycles.",
    )
    parser_test.add_argument(
        "--alpha",
        type=float,
        default=1e-4,
        help="p-value cut-off for calling a CBS breakpoint.",
    )
    parser_test.add_argument(
        "--zscore",
        type=float,
        default=5,
        help="z-score cut-off for aberration calling.",
    )
    parser_test.add_argument(
        "--beta",
        type=float,
        default=None,
        help="When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations. Beta is a number between 0 (liberal) and 1 (conservative) and is optimally close to the purity.",
    )
    parser_test.add_argument(
        "--blacklist",
        type=str,
        default=None,
        help="Blacklist that masks regions in output, structure of header-less "
        "file: chr...(/t)startpos(/t)endpos(/n)",
    )
    parser_test.add_argument(
        "--gender",
        type=str,
        choices=["F", "M"],
        help="Force WisecondorX to analyze this case as a male (M) or a female (F)",
    )
    parser_test.add_argument(
        "--ylim",
        type=str,
        default="def",
        help="y-axis limits for plotting. e.g. [-2,2]",
    )
    parser_test.add_argument(
        "--bed",
        action="store_true",
        help="Outputs tab-delimited .bed files, containing the most important information",
    )
    parser_test.add_argument("--plot", action="store_true", help="Outputs .png plots")
    parser_test.add_argument(
        "--cairo",
        action="store_true",
        help="Uses cairo bitmap type for plotting. Might be necessary for certain setups.",
    )
    parser_test.add_argument(
        "--add-plot-title",
        action="store_true",
        help="Add the output name as plot title",
    )
    parser_test.add_argument(
        "--seed", type=int, default=None, help="Seed for segmentation algorithm"
    )
    parser_test.set_defaults(func=tool_test)

    args = parser.parse_args(sys.argv[1:])

    logging.basicConfig(
        format="[%(levelname)s - %(asctime)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=getattr(logging, args.loglevel.upper(), None),
    )
    logging.debug("args are: {}".format(args))
    args.func(args)


if __name__ == "__main__":
    main()
