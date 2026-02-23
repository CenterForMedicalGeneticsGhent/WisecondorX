#!/usr/bin/env python3

import logging
import os
import sys
import warnings

import numpy as np
import typer
from typing import Optional, List, Dict, Any
from pathlib import Path

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

app = typer.Typer(
    name="wisecondorx",
    help="WisecondorX: Copy Number Aberration detection from shallow Whole Genome Sequencing data.",
    add_completion=False,
)


def setup_logging(loglevel: str = "INFO") -> None:
    logging.basicConfig(
        format="[%(levelname)s - %(asctime)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=getattr(logging, loglevel.upper(), None),
    )


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


@app.command("convert")
def tool_convert(
    infile: Path = typer.Argument(
        ..., help="aligned reads input for conversion (.bam or .cram)"
    ),
    outfile: Path = typer.Argument(..., help="Output .npz file"),
    reference: Optional[str] = typer.Option(
        None,
        "-r",
        "--reference",
        help="Fasta reference to be used during cram conversion",
    ),
    binsize: int = typer.Option(5e3, "--binsize", help="Bin size (bp)"),
    normdup: bool = typer.Option(
        False, "--normdup", help="Do not remove duplicates", is_flag=True
    ),
) -> None:
    """
    Convert and filter aligned reads to .npz format.
    """

    # argument sanity check
    # check if infile exists and has an index
    if not (infile.exists() and infile.is_file()):
        logging.error(f"Input file {infile} does not exist or is not a file.")
        sys.exit(1)
    if infile.suffix == ".bam":
        if (
            not Path(infile, ".bai").exists()
            and not Path(infile, ".csi").exists()
        ):
            logging.error(
                "Bam inputs need to have a 'bai' or 'csi' index present. Run 'samtools index {f}' to generate the index."
            )
    elif infile.suffix == ".cram":
        if not Path(infile, ".crai").exists():
            logging.error(
                "Cram inputs need to have a 'crai' index present. Run 'samtools index {f}' to generate the index."
            )
        if not reference:
            logging.error(
                "Cram inputs need a reference fasta provided through the '--reference' flag."
            )
        elif not reference.exists():
            logging.error(f"Fasta reference file {reference} does not exist.")

    logging.info("Starting conversion")

    sample, qual_info = convert_reads(
        infile=infile, binsize=binsize, normdup=normdup, reference=reference
    )
    np.savez_compressed(
        outfile, binsize=binsize, sample=sample, quality=qual_info
    )

    logging.info("Finished conversion")


@app.command("newref")
def tool_newref(
    infiles: List[Path] = typer.Argument(
        ...,
        help="Path to all reference data files (e.g. path/to/reference/*.npz)",
    ),
    outfile: Path = typer.Argument(
        ...,
        help="Path and filename for the reference output (e.g. path/to/myref.npz)",
    ),
    nipt: bool = typer.Option(
        False, "--nipt", help="Use flag for NIPT", is_flag=True
    ),
    yfrac: float = typer.Option(
        None,
        "--yfrac",
        help="Use to manually set the Y read fraction cutoff, which defines gender",
    ),
    plotyfrac: Path = typer.Option(
        None,
        "--plotyfrac",
        help="Path to yfrac .png plot for optimization; software will stop after plotting",
    ),
    refsize: int = typer.Option(
        300, "--refsize", help="Amount of reference locations per target"
    ),
    binsize: int = int(
        1e5
    ),  # Cannot use scientific notation as default for Typer Option int
    cpus: int = typer.Option(
        1, "--cpus", help="Use multiple cores to find reference bins"
    ),
) -> None:
    """
    Create a new reference using healthy reference samples.
    """
    # Fix the binsize type from typer, because 1e5 is float
    binsize_val = int(binsize)

    logging.info("Creating new reference")

    if yfrac is not None:
        if yfrac < 0 or yfrac > 1:
            logging.critical(
                "Parameter --yfrac should be a positive number lower than or equal to 1"
            )
            sys.exit(1)

    split_path = list(os.path.split(outfile))
    if split_path[-1].endswith(".npz"):
        split_path[-1] = split_path[-1][:-4]

    basepath = os.path.join(split_path[0], split_path[1])
    prepfile = f"{basepath}_prep.npz"
    prepdatafile = f"{basepath}_prep_data.npy"
    partfile = f"{basepath}_part"
    tmpoutfile_A = f"{basepath}.tmp.A.npz"
    tmpoutfile_F = f"{basepath}.tmp.F.npz"
    tmpoutfile_M = f"{basepath}.tmp.M.npz"

    samples: list[dict[str, np.ndarray]] = []
    logging.info("Importing data ...")
    for infile in infiles:
        logging.info(f"Loading: {infile}")
        npzdata = np.load(infile, encoding="latin1", allow_pickle=True)
        sample = npzdata["sample"].item()
        b_size = int(npzdata["binsize"])
        logging.info(f"Binsize: {int(b_size)}")
        samples.append(scale_sample(sample, b_size, binsize_val))

    samples_array = np.array(samples)
    genders, trained_cutoff = train_gender_model(
        samples=samples_array, yfrac=yfrac, plotyfrac=plotyfrac
    )

    if genders.count("F") < 5 and nipt:
        logging.warning(
            "A NIPT reference should have at least 5 female feti samples. Removing --nipt flag."
        )
        nipt = False

    if not nipt:
        for i, sample in enumerate(samples_array):
            samples_array[i] = gender_correct(sample, genders[i])

    total_mask, bins_per_chr = get_mask(samples_array)
    if genders.count("F") > 4:
        mask_F, _ = get_mask(samples_array[np.array(genders) == "F"])
        total_mask = total_mask & mask_F
    if genders.count("M") > 4 and not nipt:
        mask_M, _ = get_mask(samples_array[np.array(genders) == "M"])
        total_mask = total_mask & mask_M

    outfiles_list: List[str] = []

    if len(genders) > 9:
        logging.info("Starting autosomal reference creation ...")
        outfiles_list.append(tmpoutfile_A)
        tool_newref_prep(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            binsize=binsize_val,
            samples=samples_array,
            gender="A",
            mask=total_mask,
            bins_per_chr=bins_per_chr,
        )
        logging.info("This might take a while ...")
        tool_newref_main(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            partfile=partfile,
            tmpoutfile=tmpoutfile_A,
            refsize=refsize,
            cpus=cpus,
        )
    else:
        logging.critical(
            "Provide at least 10 samples to enable the generation of a reference."
        )
        sys.exit(1)

    if genders.count("F") > 4:
        logging.info("Starting female gonosomal reference creation ...")
        outfiles_list.append(tmpoutfile_F)
        tool_newref_prep(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            binsize=binsize_val,
            samples=samples_array[np.array(genders) == "F"],
            gender="F",
            mask=total_mask,
            bins_per_chr=bins_per_chr,
        )
        logging.info("This might take a while ...")
        tool_newref_main(
            prepdatafile=prepdatafile,
            prepfile=prepfile,
            partfile=partfile,
            tmpoutfile=tmpoutfile_F,
            refsize=refsize,
            cpus=1,
        )
    else:
        logging.warning(
            "Provide at least 5 female samples to enable normalization of female gonosomes."
        )

    if not nipt:
        if genders.count("M") > 4:
            logging.info("Starting male gonosomal reference creation ...")
            outfiles_list.append(tmpoutfile_M)
            tool_newref_prep(
                prepdatafile=prepdatafile,
                prepfile=prepfile,
                binsize=binsize_val,
                samples=samples_array[np.array(genders) == "M"],
                gender="M",
                mask=total_mask,
                bins_per_chr=bins_per_chr,
            )
            tool_newref_main(
                prepdatafile=prepdatafile,
                prepfile=prepfile,
                partfile=partfile,
                tmpoutfile=tmpoutfile_M,
                refsize=refsize,
                cpus=1,
            )
        else:
            logging.warning(
                "Provide at least 5 male samples to enable normalization of male gonosomes."
            )

    tool_newref_merge(
        outfile=outfile,
        nipt=nipt,
        outfiles=outfiles_list,
        trained_cutoff=trained_cutoff,
    )

    logging.info("Finished creating reference")


@app.command("gender")
def output_gender(
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


@app.command("predict")
def tool_predict(
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
        is_flag=True,
    ),
    plot: bool = typer.Option(
        False, "--plot", help="Outputs .png plots", is_flag=True
    ),
    cairo: bool = typer.Option(
        False,
        "--cairo",
        help="Uses cairo bitmap type for plotting. Might be necessary for certain setups.",
        is_flag=True,
    ),
    add_plot_title: bool = typer.Option(
        False,
        "--add-plot-title",
        help="Add the output name as plot title",
        is_flag=True,
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


def main() -> None:
    app()


if __name__ == "__main__":
    main()
