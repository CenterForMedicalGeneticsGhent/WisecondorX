import logging
import psutil
import sys
from typing import Optional

import numpy as np
import typer
from pathlib import Path
from dataclasses import dataclass
import wisecondorx_rs

"""
Converts aligned reads file to numpy array by transforming
individual reads to counts per bin.
"""


@dataclass
class BamQualityInfo:
    mapped: int
    unmapped: int
    no_coordinate: int
    filter_rmdup: int
    filter_mapq: int
    pre_retro: int
    post_retro: int
    pair_fail: int


def wcx_convert(
    infile: Path = typer.Argument(
        ..., help="aligned reads input for conversion (.bam, .cram, .sam)"
    ),
    outfile: Path = typer.Argument(..., help="Output .npz file"),
    reference: Optional[str] = typer.Option(
        None,
        "-r",
        "--reference",
        help="Fasta reference to be used during cram conversion",
    ),
    binsize: int = typer.Option(5000, "--binsize", help="Bin size (bp)"),
    rmdup: bool = typer.Option(True, "--rmdup", help="Remove duplicates"),
    threads: int = typer.Option(
        psutil.cpu_count(logical=False),
        "--threads",
        help="Number of threads for parallel processing",
    ),
) -> None:
    """
    Convert and filter aligned reads to .npz format.
    """

    if not (infile.exists() and infile.is_file()):
        logging.error(f"Input file {infile} does not exist or is not a file.")
        sys.exit(1)

    logging.info(
        "Importing and converting aligned reads with wisecondorx-rs engines ..."
    )

    sample_dict, qual_info = wisecondorx_rs.wcx_convert_core(
        str(infile), reference, binsize, rmdup, threads
    )

    logging.info("Saving processed datasets to npz ...")

    rust_qual = qual_info
    qual_info = BamQualityInfo(
        mapped=rust_qual.mapped,
        unmapped=rust_qual.unmapped,
        no_coordinate=rust_qual.no_coordinate,
        filter_rmdup=rust_qual.filter_rmdup,
        filter_mapq=rust_qual.filter_mapq,
        pre_retro=rust_qual.pre_retro,
        post_retro=rust_qual.post_retro,
        pair_fail=rust_qual.pair_fail,
    )

    # We cast to standard python dataclass fields equivalent or keep the Rust struct bound
    np.savez_compressed(
        outfile,
        binsize=binsize,
        sample=sample_dict,
        quality=qual_info,
    )

    logging.info("Finished conversion")
