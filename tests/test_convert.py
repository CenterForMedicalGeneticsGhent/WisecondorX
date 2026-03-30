import unittest
import os
import shutil
import tempfile
from pathlib import Path
import numpy as np
import pysam
from wisecondorx.convert import wcx_convert


class TestConvert(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.bam_path = os.path.join(self.test_dir, "test.bam")
        self.create_dummy_bam(self.bam_path)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def create_dummy_bam(self, path):
        header = {
            "HD": {"VN": "1.0"},
            "SQ": [
                {"LN": 10000, "SN": "1"},
                {"LN": 10000, "SN": "chr2"},
                {"LN": 10000, "SN": "X"},
                {"LN": 10000, "SN": "Y"},
            ],
        }

        with pysam.AlignmentFile(path, "wb", header=header) as outf:
            # Chromosome 1: 10 reads in first bin (0-4999), 5 in second (5000-9999)
            for i in range(10):
                a = pysam.AlignedSegment()
                a.query_name = f"read_1_{i}"
                a.query_sequence = "AGCT" * 25
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 100 + i  # Unique positions
                a.mapping_quality = 20
                a.cigar = ((0, 100),)
                outf.write(a)
            for i in range(5):
                a = pysam.AlignedSegment()
                a.query_name = f"read_1_5000_{i}"
                a.query_sequence = "AGCT" * 25
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 6000 + i  # Unique positions
                a.mapping_quality = 20
                a.cigar = ((0, 100),)
                outf.write(a)

            # Chromosome chr2: 3 reads
            for i in range(3):
                a = pysam.AlignedSegment()
                a.query_name = f"read_2_{i}"
                a.query_sequence = "AGCT" * 25
                a.flag = 0
                a.reference_id = 1
                a.reference_start = 500 + i  # Unique positions
                a.mapping_quality = 20
                a.cigar = ((0, 100),)
                outf.write(a)

            # Chromosome X (23): 2 reads
            for i in range(2):
                a = pysam.AlignedSegment()
                a.query_name = f"read_X_{i}"
                a.query_sequence = "AGCT" * 25
                a.flag = 0
                a.reference_id = 2
                a.reference_start = 500 + i  # Unique positions
                a.mapping_quality = 20
                a.cigar = ((0, 100),)
                outf.write(a)

        pysam.index(path)

    def test_wcx_convert_basic(self):
        prefix = os.path.join(self.test_dir, "output")
        wcx_convert(Path(self.bam_path), Path(prefix), binsize=5000)

        out_npz = prefix + ".npz"
        self.assertTrue(os.path.exists(out_npz))

        data = np.load(out_npz, allow_pickle=True)
        self.assertEqual(data["binsize"], 5000)

        reads_per_bin = data["reads_per_bin"].item()
        # Chr 1
        self.assertEqual(reads_per_bin["1"][0], 10)
        self.assertEqual(reads_per_bin["1"][1], 5)
        # Chr 2 (was chr2 in BAM)
        self.assertEqual(reads_per_bin["2"][0], 3)
        # Chr X (was X in BAM, should be 23)
        self.assertEqual(reads_per_bin["23"][0], 2)

        quality = data["quality"].item()
        self.assertEqual(quality["post_retro"], 20)  # 10 + 5 + 3 + 2

    def test_wcx_convert_rmdup(self):
        # Create a BAM with duplicates
        dup_bam = os.path.join(self.test_dir, "dup.bam")
        header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 10000, "SN": "1"}]}
        with pysam.AlignmentFile(dup_bam, "wb", header=header) as outf:
            for i in range(2):
                a = pysam.AlignedSegment()
                a.query_name = "read_dup"
                a.query_sequence = "AGCT" * 25
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 100  # Same position
                a.mapping_quality = 20
                a.cigar = ((0, 100),)
                outf.write(a)
        pysam.index(dup_bam)

        prefix = os.path.join(self.test_dir, "output_dup")
        # normdup=False (default) should remove duplicates
        wcx_convert(Path(dup_bam), Path(prefix), binsize=5000, normdup=False)
        data = np.load(prefix + ".npz", allow_pickle=True)
        reads_per_bin = data["reads_per_bin"].item()
        self.assertEqual(reads_per_bin["1"][0], 1)
        self.assertEqual(data["quality"].item()["filter_rmdup"], 1)

        # normdup=True should keep duplicates
        wcx_convert(Path(dup_bam), Path(prefix), binsize=5000, normdup=True)
        data = np.load(prefix + ".npz", allow_pickle=True)
        reads_per_bin = data["reads_per_bin"].item()
        self.assertEqual(reads_per_bin["1"][0], 2)
        self.assertEqual(data["quality"].item()["filter_rmdup"], 0)

    def test_wcx_convert_missing_file(self):
        with self.assertRaises(SystemExit) as cm:
            wcx_convert(Path("non_existent.bam"), Path("out"))
        self.assertEqual(cm.exception.code, 1)

    def test_wcx_convert_cram_no_ref(self):
        # Use a real BAM but name it .cram to trigger my logic check without pysam failing on headers
        cram_path = os.path.join(self.test_dir, "test_fake.cram")
        shutil.copy(self.bam_path, cram_path)
        with self.assertRaises(SystemExit) as cm:
            wcx_convert(Path(cram_path), Path("out"))
        self.assertEqual(cm.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
