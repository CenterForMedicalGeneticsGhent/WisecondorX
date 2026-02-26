import json
from pathlib import Path
from tempfile import TemporaryDirectory
from wisecondorx.predict import _get_processed_cbs


def test_predict_stub():
    """Stub test for predict module."""
    import wisecondorx.predict

    assert wisecondorx.predict is not None


def test_get_processed_cbs():
    with TemporaryDirectory() as tmpdir:
        cbs_json = Path(tmpdir) / "cbs_output.json"

        # Create some dummy CBS output data
        dummy_data = [
            {"chromosome": 1, "start": 1000, "end": 2000, "ratio": 1.5},
            {"chromosome": "2", "start": "3000", "end": "4000", "ratio": -0.5},
            {"chromosome": 22, "start": 5000, "end": 6000, "ratio": 0.0},
        ]

        with open(cbs_json, "w") as f:
            json.dump(dummy_data, f)

        # Call the function
        segments = _get_processed_cbs(cbs_json)

        # Verify the output
        assert len(segments) == 3
        # chromosome is 0-indexed in the output
        assert segments[0] == (0, 1000, 2000, 1.5)
        assert segments[1] == (1, 3000, 4000, -0.5)
        assert segments[2] == (21, 5000, 6000, 0.0)


def test_get_processed_cbs_empty():
    with TemporaryDirectory() as tmpdir:
        cbs_json = Path(tmpdir) / "cbs_empty.json"

        with open(cbs_json, "w") as f:
            json.dump([], f)

        segments = _get_processed_cbs(cbs_json)
        assert len(segments) == 0
