# WisecondorX

import logging

import sys

import numpy as np
import pysam
import typer
from pathlib import Path


from typing import Annotated


def wcx_convert(
    infile: Annotated[
        Path,
        typer.Argument(
            help="aligned reads input for conversion (.bam or .cram)",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ],
    prefix: Annotated[Path, typer.Argument(help="Output prefix")],
    reference: Annotated[
        Path,
        typer.Option(
            "-r",
            "--reference",
            help="Fasta reference to be used during cram conversion",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ] = None,
    binsize: Annotated[
        int, typer.Option("--binsize", help="Bin size (bp)")
    ] = 5000,
    normdup: Annotated[
        bool, typer.Option("--normdup", help="Do not remove duplicates")
    ] = False,
) -> None:
    """
    Convert and filter aligned reads to .npz format.
    """

    # check if infile exists and has an index
    if not (infile.exists() and infile.is_file()):
        logging.error(f"Input file {infile} does not exist or is not a file.")
        sys.exit(1)

    logging.info("Importing data ...")

    reads_per_chromosome_bin: dict[str, np.ndarray] = {
        str(i): None for i in range(1, 25)
    }

    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0

    logging.info("Converting aligned reads")

    mode = "rb" if infile.suffix == ".bam" else "rc"
    ref_filename = str(reference) if reference else None

    if mode == "rc" and not ref_filename:
        logging.error(
            "Cram inputs need a reference fasta provided through the '--reference' flag."
        )
        sys.exit(1)

    with pysam.AlignmentFile(
        infile, mode, reference_filename=ref_filename
    ) as reads_file:
        for index, chr in enumerate(reads_file.references):
            larp = -1
            larp2 = -1
            chr_name = chr
            if chr_name[:3].lower() == "chr":
                chr_name = chr_name[3:]

            if chr_name == "X":
                chr_name = "23"
            elif chr_name == "Y":
                chr_name = "24"

            if chr_name not in reads_per_chromosome_bin:
                continue

            logging.info(
                f"Working at {chr}; processing {int(reads_file.lengths[index] / float(binsize) + 1)} bins"
            )
            counts = np.zeros(
                int(reads_file.lengths[index] / float(binsize) + 1),
                dtype=np.int32,
            )
            bam_chr = reads_file.fetch(chr)

            for read in bam_chr:
                if read.is_paired:
                    if not read.is_proper_pair:
                        reads_pairf += 1
                        continue
                    if (
                        not normdup
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
                    if not normdup and larp == read.pos:
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

        qual_info: dict[str, int] = {
            "mapped": reads_file.mapped,
            "unmapped": reads_file.unmapped,
            "no_coordinate": reads_file.nocoordinate,
            "filter_rmdup": reads_rmdup,
            "filter_mapq": reads_mapq,
            "pre_retro": reads_seen,
            "post_retro": reads_kept,
            "pair_fail": reads_pairf,
        }

    np.savez_compressed(
        Path(f"{prefix}.npz"),
        binsize=binsize,
        reads_per_bin=reads_per_chromosome_bin,
        quality=qual_info,
    )

    logging.info("Finished conversion")
