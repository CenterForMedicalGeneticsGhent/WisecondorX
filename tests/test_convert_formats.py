import json
import tempfile
import unittest
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


from wisecondorx.convert import (
    ConvertOutput,
    load_convert_output,
    wcx_convert,
)


class FakeRead:
    """Mock read for pysam."""

    def __init__(
        self,
        reference_start: int,
        mapping_quality: int = 10,
        is_paired: bool = False,
        is_proper_pair: bool = True,
        next_reference_start: int = 0,
    ):
        self.reference_start = reference_start
        self.mapping_quality = mapping_quality
        self.is_paired = is_paired
        self.is_proper_pair = is_proper_pair
        self.next_reference_start = next_reference_start


class FakeAlignmentFile:
    """Mock AlignmentFile for pysam."""

    def __init__(self, *_args: Any, **_kwargs: Any):
        self.references = ["chr1", "chrX", "chrM"]
        self.lengths = [300, 200, 100]
        self.mapped = 3
        self.unmapped = 1
        self.nocoordinate = 0

    def fetch(self, chromosome: str) -> Any:
        reads = {
            "chr1": [
                FakeRead(reference_start=0, mapping_quality=10),
                FakeRead(reference_start=100, mapping_quality=0),
            ],
            "chrX": [FakeRead(reference_start=0, mapping_quality=10)],
        }
        return iter(reads.get(chromosome, []))

    def close(self) -> None:
        pass

    def __enter__(self) -> "FakeAlignmentFile":
        return self

    def __exit__(self, *_exc: Any) -> None:
        pass


class TestConvertOutputFormats(unittest.TestCase):
    """Tests for coordinate/format conversion in WisecondorX convert."""

    def test_convert_writes_legacy_npz_contract(self):
        """Verify wcx_convert output conforms to the legacy .npz contract."""
        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / "in.bam"
            infile.touch()
            prefix = Path(tmpdir) / "sample"

            with (
                patch(
                    "wisecondorx.convert.pysam.AlignmentFile",
                    FakeAlignmentFile,
                ),
                patch("wisecondorx.convert.pysam.index"),
            ):
                wcx_convert(
                    infile=infile,
                    prefix=prefix,
                    binsize=100,
                    normdup=True,
                    threads=1,
                    out_format=ConvertOutput.NPZ,
                )

            npz_path = Path(f"{prefix}.npz")
            self.assertTrue(npz_path.exists())

            with np.load(npz_path, allow_pickle=True) as out:
                self.assertIn("binsize", out)
                self.assertIn("sample", out)
                self.assertIn("quality", out)

                sample = out["sample"].item()
                quality = out["quality"].item()

                self.assertEqual(set(sample.keys()), {"1", "23"})
                self.assertEqual(int(sample["1"][0]), 1)
                self.assertEqual(int(sample["23"][0]), 1)
                self.assertEqual(int(quality["filter_mapq"]), 1)

    def test_convert_both_writes_npz_and_parquet(self):
        """Verify wcx_convert writes both correct .npz and .parquet files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / "in.bam"
            infile.touch()
            prefix = Path(tmpdir) / "sample"

            with (
                patch(
                    "wisecondorx.convert.pysam.AlignmentFile",
                    FakeAlignmentFile,
                ),
                patch("wisecondorx.convert.pysam.index"),
            ):
                wcx_convert(
                    infile=infile,
                    prefix=prefix,
                    binsize=100,
                    normdup=True,
                    threads=1,
                    out_format=ConvertOutput.BOTH,
                )

            self.assertTrue(Path(f"{prefix}.npz").exists())
            pq_path = Path(f"{prefix}.parquet")
            self.assertTrue(pq_path.exists())

            # Verify the written parquet metadata
            table = pq.read_table(pq_path)
            metadata = table.schema.metadata
            self.assertIsNotNone(metadata)
            self.assertIn(b"wisecondorx.schema", metadata)
            self.assertIn(b"wisecondorx.binsize", metadata)
            self.assertIn(b"wisecondorx.quality", metadata)
            self.assertIn(b"wisecondorx.chromosomes", metadata)
            self.assertEqual(metadata[b"wisecondorx.binsize"], b"100")

            quality = json.loads(metadata[b"wisecondorx.quality"].decode())
            self.assertEqual(int(quality["filter_mapq"]), 1)

    def test_convert_parquet_only_skips_npz(self):
        """Verify wcx_convert can output only parquet, skipping .npz."""
        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / "in.bam"
            infile.touch()
            prefix = Path(tmpdir) / "sample"

            with (
                patch(
                    "wisecondorx.convert.pysam.AlignmentFile",
                    FakeAlignmentFile,
                ),
                patch("wisecondorx.convert.pysam.index"),
            ):
                wcx_convert(
                    infile=infile,
                    prefix=prefix,
                    binsize=100,
                    normdup=True,
                    threads=1,
                    out_format=ConvertOutput.PARQUET,
                )

            self.assertFalse(Path(f"{prefix}.npz").exists())
            self.assertTrue(Path(f"{prefix}.parquet").exists())

    def test_convert_invalid_out_format_raises(self):
        """Verify an invalid output format type/value raises an OSError."""
        with self.assertRaises(OSError):
            wcx_convert(
                infile=Path("missing.bam"),
                prefix=Path("sample"),
                binsize=100,
                normdup=True,
                threads=1,
                out_format="invalid",  # type: ignore
            )

    def test_load_convert_output_reads_npz_and_parquet(self):
        """Verify load_convert_output loads both .npz and .parquet formats."""
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "sample"

            # 1. Test .npz loading
            np.savez_compressed(
                Path(f"{prefix}.npz"),
                binsize=100,
                sample=np.array(
                    {"1": np.array([1, 2]), "23": np.array([3])}, dtype=object
                ),
                quality=np.array({"foo": 1}, dtype=object),
            )

            npz_sample, npz_binsize = load_convert_output(
                Path(f"{prefix}.npz")
            )
            self.assertEqual(npz_binsize, 100)
            self.assertEqual(set(npz_sample.keys()), {"1", "23"})

            # 2. Test .parquet loading (write real parquet using pyarrow)
            counts_df = pd.DataFrame(
                {"chr": [1, 1, 23], "bin": [0, 1, 0], "count": [1, 0, 1]}
            )
            table = pa.Table.from_pandas(counts_df, preserve_index=False)
            metadata = {
                b"wisecondorx.binsize": b"100",
                b"wisecondorx.quality": b"{}",
                b"wisecondorx.chromosomes": b"[]",
            }
            table = table.replace_schema_metadata(metadata)

            pq_path = Path(f"{prefix}.parquet")
            pq.write_table(table, pq_path)

            parquet_sample, parquet_binsize = load_convert_output(pq_path)
            self.assertEqual(parquet_binsize, 100)
            self.assertEqual(set(parquet_sample.keys()), {"1", "23"})
            self.assertTrue(
                np.array_equal(parquet_sample["1"], np.array([1, 0]))
            )
            self.assertTrue(
                np.array_equal(parquet_sample["23"], np.array([1]))
            )


if __name__ == "__main__":
    unittest.main()
