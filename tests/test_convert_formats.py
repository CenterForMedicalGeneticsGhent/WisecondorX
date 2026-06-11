import tempfile
import unittest
import os
import sys
import types
import json
from typing import Any
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
)

if "pysam" not in sys.modules:
    pysam_stub: Any = types.ModuleType("pysam")
    pysam_stub.AlignmentFile = None
    pysam_stub.index = lambda *_args, **_kwargs: None
    sys.modules["pysam"] = pysam_stub

from wisecondorx.convert import wcx_convert, ConvertOutput
import wisecondorx.convert as convert_module


class _FakeRead:
    def __init__(
        self,
        reference_start,
        mapping_quality=10,
        is_paired=False,
        is_proper_pair=True,
        next_reference_start=0,
    ):
        self.reference_start = reference_start
        self.mapping_quality = mapping_quality
        self.is_paired = is_paired
        self.is_proper_pair = is_proper_pair
        self.next_reference_start = next_reference_start


class _FakeAlignmentFile:
    def __init__(self, *_args, **_kwargs):
        self.references = ["chr1", "chrX", "chrM"]
        self.lengths = [300, 200, 100]
        self.mapped = 3
        self.unmapped = 1
        self.nocoordinate = 0

    def fetch(self, chromosome):
        reads = {
            "chr1": [
                _FakeRead(reference_start=0, mapping_quality=10),
                _FakeRead(reference_start=100, mapping_quality=0),
            ],
            "chrX": [_FakeRead(reference_start=0, mapping_quality=10)],
        }
        return iter(reads.get(chromosome, []))

    def close(self):
        return None

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        return None


class _FakeArrowSchema:
    def __init__(self):
        self.metadata: Any = None


class _FakeArrowTable:
    def __init__(self):
        self.schema = _FakeArrowSchema()
        self._frame = None

    def replace_schema_metadata(self, metadata):
        self.schema.metadata = metadata
        return self

    def to_pandas(self):
        return self._frame


class _FakeArrowTableFactory:
    @staticmethod
    def from_pandas(_df, preserve_index=False):
        if preserve_index:
            raise AssertionError("preserve_index should be False")
        table = _FakeArrowTable()
        table._frame = _df
        return table


def _fake_pyarrow_modules(calls):
    pyarrow_mod: Any = types.ModuleType("pyarrow")
    pyarrow_mod.Table = _FakeArrowTableFactory
    pyarrow_mod.__path__ = []

    parquet_mod: Any = types.ModuleType("pyarrow.parquet")

    def _write_table(table, path, compression="zstd"):
        calls.append(
            {
                "path": str(path),
                "compression": compression,
                "metadata": table.schema.metadata,
            }
        )

    parquet_mod.write_table = _write_table

    def _read_table(path):
        table = _FakeArrowTable()
        table.schema.metadata = {
            b"wisecondorx.binsize": b"100",
            b"wisecondorx.quality": b"{}",
            b"wisecondorx.chromosomes": b"[]",
        }
        table._frame = __import__("pandas").DataFrame(
            {"chr": [1, 1, 23], "bin": [0, 1, 0], "count": [1, 0, 1]}
        )
        return table

    parquet_mod.read_table = _read_table
    pyarrow_mod.parquet = parquet_mod
    return {"pyarrow": pyarrow_mod, "pyarrow.parquet": parquet_mod}


def _fake_pyarrow_objects(calls):
    pyarrow_mod: Any = types.SimpleNamespace(Table=_FakeArrowTableFactory)
    parquet_mod: Any = types.SimpleNamespace()

    def _write_table(table, path, compression="zstd"):
        calls.append(
            {
                "path": str(path),
                "compression": compression,
                "metadata": table.schema.metadata,
            }
        )

    parquet_mod.write_table = _write_table
    return pyarrow_mod, parquet_mod


class TestConvertOutputFormats(unittest.TestCase):
    def test_convert_writes_legacy_npz_contract(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / "in.bam"
            infile.touch()
            prefix = Path(tmpdir) / "sample"

            with (
                patch(
                    "wisecondorx.convert.pysam.AlignmentFile",
                    _FakeAlignmentFile,
                ),
                patch("wisecondorx.convert.pysam.index") as _index,
            ):
                wcx_convert(
                    infile=infile,
                    prefix=prefix,
                    binsize=100,
                    normdup=True,
                    threads=1,
                    out_format=ConvertOutput.NPZ,
                )

            out = np.load(Path(f"{prefix}.npz"), allow_pickle=True)
            try:
                self.assertIn("binsize", out)
                self.assertIn("sample", out)
                self.assertIn("quality", out)

                sample = out["sample"].item()
                quality = out["quality"].item()

                self.assertEqual(set(sample.keys()), {"1", "23"})
                self.assertEqual(int(sample["1"][0]), 1)
                self.assertEqual(int(sample["23"][0]), 1)
                self.assertEqual(int(quality["filter_mapq"]), 1)
            finally:
                out.close()

    def test_convert_both_writes_npz_and_parquet(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / "in.bam"
            infile.touch()
            prefix = Path(tmpdir) / "sample"
            pq_calls = []
            fake_pa, fake_pq = _fake_pyarrow_objects(pq_calls)

            with (
                patch(
                    "wisecondorx.convert.pysam.AlignmentFile",
                    _FakeAlignmentFile,
                ),
                patch("wisecondorx.convert.pysam.index") as _index,
                patch.object(convert_module, "pa", fake_pa),
                patch.object(convert_module, "pq", fake_pq),
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
            self.assertEqual(len(pq_calls), 1)
            self.assertEqual(pq_calls[0]["path"], f"{prefix}.parquet")
            self.assertEqual(pq_calls[0]["compression"], "zstd")

            metadata = pq_calls[0]["metadata"]
            self.assertIn(b"wisecondorx.schema", metadata)
            self.assertIn(b"wisecondorx.binsize", metadata)
            self.assertIn(b"wisecondorx.quality", metadata)
            self.assertIn(b"wisecondorx.chromosomes", metadata)
            self.assertEqual(metadata[b"wisecondorx.binsize"], b"100")

            quality = json.loads(metadata[b"wisecondorx.quality"].decode())
            self.assertEqual(int(quality["filter_mapq"]), 1)

    def test_convert_parquet_only_skips_npz(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            infile = Path(tmpdir) / "in.bam"
            infile.touch()
            prefix = Path(tmpdir) / "sample"
            pq_calls = []
            fake_pa, fake_pq = _fake_pyarrow_objects(pq_calls)

            with (
                patch(
                    "wisecondorx.convert.pysam.AlignmentFile",
                    _FakeAlignmentFile,
                ),
                patch("wisecondorx.convert.pysam.index") as _index,
                patch.object(convert_module, "pa", fake_pa),
                patch.object(convert_module, "pq", fake_pq),
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
            self.assertEqual(len(pq_calls), 1)
            self.assertEqual(pq_calls[0]["path"], f"{prefix}.parquet")

    def test_convert_invalid_out_format_exits(self):
        with self.assertRaises(TypeError):
            wcx_convert(
                infile=Path("missing.bam"),
                prefix=Path("sample"),
                out_format=ConvertOutput.BOTH,
            )

    def test_load_convert_output_reads_npz_and_parquet(self):
        from wisecondorx.convert import load_convert_output

        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "sample"
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

            pq_calls = []
            with patch.dict(sys.modules, _fake_pyarrow_modules(pq_calls)):
                parquet_sample, parquet_binsize = load_convert_output(
                    Path(f"{prefix}.parquet")
                )

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
