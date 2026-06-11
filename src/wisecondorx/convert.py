# WisecondorX

import logging
import os
import re
import json

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
import pandas as pd
import pysam
import typer
from pathlib import Path
from enum import Enum


def load_convert_output(infile: Path) -> tuple[dict[str, np.ndarray], int]:
    """
    Load a conversion output file produced by wcx_convert.

    Supports the legacy .npz contract and the single-file parquet format.
    """

    if infile.suffix == ".npz":
        npzdata = np.load(infile, encoding="latin1", allow_pickle=True)
        try:
            sample = npzdata["sample"].item()
            binsize = int(npzdata["binsize"])
        finally:
            npzdata.close()
        return sample, binsize

    if infile.suffix == ".parquet":
        try:
            import pyarrow.parquet as pq  # type: ignore[import-not-found]
        except (ImportError, ModuleNotFoundError) as error:
            raise ValueError(
                "Parquet input requires pyarrow to be installed"
            ) from error

        table = pq.read_table(infile)
        metadata = table.schema.metadata or {}

        try:
            binsize = int(metadata[b"wisecondorx.binsize"].decode("utf-8"))
        except KeyError as error:
            raise ValueError(
                "Parquet input is missing required metadata key: binsize"
            ) from error

        df = table.to_pandas()
        sample: dict[str, np.ndarray] = {}
        for chr_num in sorted(df["chr"].drop_duplicates().astype(int)):
            chr_df = df[df["chr"] == chr_num].sort_values("bin")
            sample[str(int(chr_num))] = chr_df["count"].to_numpy()

        return sample, binsize

    raise ValueError(
        "Unsupported conversion output format: {}".format(infile.suffix)
    )


class ConvertOutput(Enum):
    BOTH = "both"
    NPZ = "npz"
    PARQUET = "parquet"


def wcx_convert(
    infile: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="aligned reads input for conversion (.bam or .cram)",
    ),
    prefix: Path = typer.Argument(..., help="Output prefix"),
    reference: Path = typer.Option(
        None,
        "-r",
        "--reference",
        help="Fasta reference to be used during cram conversion",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
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
    out_format: ConvertOutput = typer.Option(
        ConvertOutput.NPZ,
        "--out-format",
        "-o",
        case_sensitive=False,
        show_default=True,
        help='Output format, options are "npz", "parquet", or "both"',
    ),
) -> None:
    """
    Convert and filter aligned reads to .npz or parquet format.
    """

    reads_file: pysam.AlignmentFile
    # check if infile has an index
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
            raise typer.BadParameter(
                "Cram inputs need a reference fasta provided through the '--reference' flag."
            )
        if not infile.with_suffix(infile.suffix + ".crai").exists():
            logging.warning(
                f"Input file {infile} does not have an index (.crai). Indexing prior to analysis..."
            )
            pysam.index(str(infile))
        reads_file = pysam.AlignmentFile(
            str(infile), "rc", reference_filename=str(reference)
        )
    else:
        raise typer.BadParameter(
            f"Input file {infile} should have extension .bam or .cram"
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
            r"^(?:chr)?(?:[1-9]|1[1-9]|2[0-2]|[XY])$",
            chromosome,
            re.IGNORECASE,
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

    if out_format in {ConvertOutput.NPZ, ConvertOutput.BOTH}:
        np.savez_compressed(
            Path(f"{prefix}.npz"),
            binsize=binsize,
            sample=np.array(reads_per_chromosome_bin, dtype=object),
            quality=np.array(qual_info, dtype=object),
        )

    if out_format in {ConvertOutput.PARQUET, ConvertOutput.BOTH}:
        counts_frames = []
        chr_meta = []
        for chr_num in sorted(reads_per_chromosome_bin.keys(), key=int):
            counts = reads_per_chromosome_bin[chr_num]
            if counts is None:
                continue
            chr_int = int(chr_num)
            counts_frames.append(
                pd.DataFrame(
                    {
                        "chr": np.full(len(counts), chr_int, dtype=np.int16),
                        "bin": np.arange(len(counts), dtype=np.int32),
                        "count": counts.astype(np.int32, copy=False),
                    }
                )
            )
            chr_meta.append(
                {
                    "chr": chr_int,
                    "n_bins": int(len(counts)),
                    "binsize": int(binsize),
                }
            )

        counts_df = (
            pd.concat(counts_frames, ignore_index=True)
            if counts_frames
            else pd.DataFrame(columns=["chr", "bin", "count"])
        )

        table = pa.Table.from_pandas(counts_df, preserve_index=False)
        metadata = dict(table.schema.metadata or {})
        metadata.update(
            {
                b"wisecondorx.schema": b"convert-counts-v1",
                b"wisecondorx.binsize": str(int(binsize)).encode("utf-8"),
                b"wisecondorx.quality": json.dumps(
                    qual_info, separators=(",", ":")
                ).encode("utf-8"),
                b"wisecondorx.chromosomes": json.dumps(
                    chr_meta, separators=(",", ":")
                ).encode("utf-8"),
            }
        )
        table = table.replace_schema_metadata(metadata)
        pq.write_table(
            table,
            Path(f"{prefix}.parquet"),
            compression="zstd",
        )

    reads_file.close()
    logging.info("Finished conversion")
