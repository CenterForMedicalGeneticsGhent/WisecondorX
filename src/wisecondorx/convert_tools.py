# WisecondorX

import logging
import os

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import numpy as np
import re
import pysam
import typer
from pathlib import Path


def wcx_convert(
    infile: Path = typer.Argument(
        ..., help="aligned reads input for conversion (.bam or .cram)"
    ),
    prefix: Path = typer.Argument(..., help="Output prefix"),
    reference: Path = typer.Option(
        None,
        "-r",
        "--reference",
        help="Fasta reference to be used during cram conversion",
    ),
    binsize: int = typer.Option(5000, "--binsize", help="Bin size (bp)"),
    normdup: bool = typer.Option(
        False, "--normdup", help="Do not remove duplicates"
    ),
    threads: int = typer.Option(
        os.cpu_count() or 1,
        "--threads",
        "-t",
        min=1,
        help="Number of threads for per-chromosome parallel conversion",
    ),
) -> None:
    """
    Convert and filter aligned reads to .npz format.
    """

    reads_file: pysam.AlignmentFile
    # check if infile exists and has an index
    if not (infile.exists() and infile.is_file()):
        logging.error(f"Input file {infile} does not exist or is not a file.")
        sys.exit(1)
    if infile.suffix == ".bam":
        if (
            not infile.with_suffix(infile.suffix + ".bai").exists()
            and not infile.with_suffix(infile.suffix + ".csi").exists()
        ):
            logging.warning(
                f"Input file {infile} does not have an index (.bai, .csi). Indexing prior to analysis..."
            )
            pysam.index(str(infile))
        reads_file = pysam.AlignmentFile(str(infile), "rb")
    elif infile.suffix == ".cram":
        if not reference:
            logging.error(
                "Cram inputs need a reference fasta provided through the '--reference' flag."
            )
        elif not reference.exists():
            logging.fatal(f"Fasta reference file {reference} does not exist.")
            sys.exit(1)
        if not infile.with_suffix(infile.suffix + ".crai").exists():
            logging.warning(
                f"Input file {infile} does not have an index (.crai). Indexing prior to analysis..."
            )
            pysam.index(str(infile))
        reads_file = pysam.AlignmentFile(
            str(infile), "rc", reference_filename=str(reference)
        )

    logging.info("Importing data ...")

    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0

    logging.info("Converting aligned reads")

    mode = "rb" if infile.suffix == ".bam" else "rc"
    ref_filename: Optional[str] = str(reference) if reference else None

    def process_chromosome(
        chromosome: str, chromosome_length: int, output_key: str
    ) -> tuple[str, np.ndarray, int, int, int, int, int]:
        local_reads_seen = 0
        local_reads_mapq = 0
        local_reads_rmdup = 0
        local_reads_pairf = 0
        larp = -1
        larp2 = -1

        logging.info(
            "Working at {}; processing {} bins".format(
                chromosome, int(chromosome_length / float(binsize) + 1)
            )
        )
        counts = np.zeros(
            int(chromosome_length / float(binsize) + 1),
            dtype=np.int32,
        )

        alignment_kwargs = {}
        if ref_filename:
            alignment_kwargs["reference_filename"] = ref_filename

        with pysam.AlignmentFile(
            str(infile), mode, **alignment_kwargs
        ) as local_reads_file:
            bam_chr = local_reads_file.fetch(chromosome)
            for read in bam_chr:
                read_start = read.reference_start
                if read_start < 0:
                    continue

                if read.is_paired:
                    if not read.is_proper_pair:
                        local_reads_pairf += 1
                        continue
                    if (
                        not normdup
                        and larp == read_start
                        and larp2 == read.next_reference_start
                    ):
                        local_reads_rmdup += 1
                    else:
                        if read.mapping_quality >= 1:
                            location = read_start / binsize
                            counts[int(location)] += 1
                        else:
                            local_reads_mapq += 1

                    larp2 = read.next_reference_start
                    local_reads_seen += 1
                    larp = read_start
                else:
                    if not normdup and larp == read_start:
                        local_reads_rmdup += 1
                    else:
                        if read.mapping_quality >= 1:
                            location = read_start / binsize
                            counts[int(location)] += 1
                        else:
                            local_reads_mapq += 1

                    local_reads_seen += 1
                    larp = read_start

        local_reads_kept = int(sum(counts))
        return (
            output_key,
            counts,
            local_reads_seen,
            local_reads_kept,
            local_reads_mapq,
            local_reads_rmdup,
            local_reads_pairf,
        )

    jobs: list[tuple[str, int, str]] = []
    reads_per_chromosome_bin: dict[str, Optional[np.ndarray]] = dict()
    for index, chromosome in enumerate(reads_file.references):
        if re.match(
            r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y?)$", chromosome, re.IGNORECASE
        ):
            chr_num = chromosome.lower().lstrip("chr")
            if chr_num == "x":
                chr_num = "23"
            elif chr_num == "y":
                chr_num = "24"
            reads_per_chromosome_bin[chr_num] = None
            jobs.append((chromosome, reads_file.lengths[index], chr_num))

    workers = min(len(jobs), threads)
    if workers > 0:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(process_chromosome, chrom, chrom_len, out_key)
                for chrom, chrom_len, out_key in jobs
            ]

            for future in as_completed(futures):
                (
                    chr_num,
                    counts,
                    local_reads_seen,
                    local_reads_kept,
                    local_reads_mapq,
                    local_reads_rmdup,
                    local_reads_pairf,
                ) = future.result()

                reads_per_chromosome_bin[chr_num] = counts
                reads_seen += local_reads_seen
                reads_kept += local_reads_kept
                reads_mapq += local_reads_mapq
                reads_rmdup += local_reads_rmdup
                reads_pairf += local_reads_pairf

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
        sample=reads_per_chromosome_bin,
        quality=qual_info,
    )

    logging.info("Finished conversion")
