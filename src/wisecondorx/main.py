#!/usr/bin/env python3

import logging
import os
import sys
import warnings
from typing import Annotated
import numpy as np
import typer
import argparse
from pathlib import Path
from wisecondorx.convert_tools import wcx_convert
from wisecondorx.predict_tools import wcx_gender, wcx_predict
from wisecondorx.newref_control import (
    tool_newref_prep,
    tool_newref_main,
    tool_newref_merge,
)
from wisecondorx.newref_tools import train_gender_model, get_mask
from wisecondorx.overall_tools import gender_correct, scale_sample
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
