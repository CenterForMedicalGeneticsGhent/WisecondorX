# WisecondorX

import logging

import sys
from typing import Optional

import numpy as np
import pysam
import typer
from pathlib import Path
from dataclasses import dataclass

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
        ..., help="aligned reads input for conversion (.bam or .cram)"
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
) -> None:
    """
    Convert and filter aligned reads to .npz format.
    """

    reads_file: pysam.AlignmentFile = None
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
        reads_file = pysam.AlignmentFile(infile, "rb")
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
        reads_file = pysam.AlignmentFile(
            infile, "rc", reference_filename=reference
        )

    logging.info("Importing data ...")

    reads_per_chromosome_bin: dict[str, np.ndarray] = dict()
    for chr in range(1, 25):
        reads_per_chromosome_bin[str(chr)] = None

    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0
    larp = -1
    larp2 = -1

    logging.info("Converting aligned reads ... This might take a while ...")

    for index, chr in enumerate(reads_file.references):
        chr_name = chr
        if chr_name[:3].lower() == "chr":
            chr_name = chr_name[3:]
        if (
            chr_name not in reads_per_chromosome_bin
            and chr_name != "X"
            and chr_name != "Y"
        ):
            continue

        logging.info(
            "Working at {}; processing {} bins".format(
                chr, int(reads_file.lengths[index] / float(binsize) + 1)
            )
        )
        counts = np.zeros(
            int(reads_file.lengths[index] / float(binsize) + 1),
            dtype=np.int32,
        )
        bam_chr = reads_file.fetch(chr)

        if chr_name == "X":
            chr_name = "23"
        if chr_name == "Y":
            chr_name = "24"

        for read in bam_chr:
            if read.is_paired:
                if not read.is_proper_pair:
                    reads_pairf += 1
                    continue
                if (
                    rmdup
                    and larp == read.pos
                    and larp2 == read.next_reference_start
                ):
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                larp2 = read.next_reference_start
                reads_seen += 1
                larp = read.pos
            else:
                if rmdup and larp == read.pos:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= 1:
                        location = read.pos / binsize
                        counts[int(location)] += 1
                    else:
                        reads_mapq += 1

                reads_seen += 1
                larp = read.pos

        reads_per_chromosome_bin[chr_name] = counts
        reads_kept += sum(counts)

    qual_info = BamQualityInfo(
        mapped=reads_file.mapped,
        unmapped=reads_file.unmapped,
        no_coordinate=reads_file.nocoordinate,
        filter_rmdup=reads_rmdup,
        filter_mapq=reads_mapq,
        pre_retro=reads_seen,
        post_retro=reads_kept,
        pair_fail=reads_pairf,
    )

    np.savez_compressed(
        outfile,
        binsize=binsize,
        sample=reads_per_chromosome_bin,
        quality=qual_info,
    )

    logging.info("Finished conversion")
