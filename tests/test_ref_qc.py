import unittest

from wisecondorx.ref_qc import _verdict_m


class TestRefQc(unittest.TestCase):
    def test_ref_qc_male_verdict_warns_on_low_chry_usable_pct(self):
        metrics = {
            "n_valid": 1000,
            "n_low_refs": 0,
            "mean_of_means": 1.0,
            "outlier_pct": 0.0,
            "chrY": {
                "n_valid": 100,
                "mean_of_means": 1.0,
                "usable_pct": 40.0,
            },
        }
        verdict, msg = _verdict_m(metrics)
        self.assertEqual(verdict, "WARN")
        self.assertIn("usable chrY bins", msg)

    def test_ref_qc_male_verdict_fails_on_very_low_chry_usable_pct(self):
        metrics = {
            "n_valid": 1000,
            "n_low_refs": 0,
            "mean_of_means": 1.0,
            "outlier_pct": 0.0,
            "chrY": {
                "n_valid": 100,
                "mean_of_means": 1.0,
                "usable_pct": 10.0,
            },
        }
        verdict, msg = _verdict_m(metrics)
        self.assertEqual(verdict, "FAIL")
        self.assertIn("usable chrY bins", msg)


if __name__ == "__main__":
    unittest.main()

