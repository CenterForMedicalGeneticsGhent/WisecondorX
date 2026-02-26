import numpy as np
import pytest
from unittest.mock import patch

from wisecondorx.utils import Sex, scale_sample, sex_correct


def test_scale_sample_same_binsize():
    sample = {"1": np.array([1, 2, 3])}
    result = scale_sample(sample, 1000, 1000)
    assert result is sample  # Should return the exact same object
    np.testing.assert_array_equal(result["1"], sample["1"])


def test_scale_sample_invalid_binsize():
    sample = {"1": np.array([1, 2, 3])}

    # Target smaller than source
    with patch("wisecondorx.utils.logging.critical") as mock_logging:
        with pytest.raises(SystemExit) as exc_info:
            scale_sample(sample, 1000, 500)
        assert exc_info.value.code == 1
        mock_logging.assert_called_once()

    # Source zero
    with patch("wisecondorx.utils.logging.critical"):
        with pytest.raises(SystemExit):
            scale_sample(sample, 0, 1000)

    # Target zero
    with patch("wisecondorx.utils.logging.critical"):
        with pytest.raises(SystemExit):
            scale_sample(sample, 1000, 0)

    # Target not a multiple of source
    with patch("wisecondorx.utils.logging.critical"):
        with pytest.raises(SystemExit):
            scale_sample(sample, 1000, 1500)


def test_scale_sample_valid():
    sample = {"1": np.array([1, 2, 3, 4, 5]), "2": np.array([10, 20])}
    # scale is 2
    result = scale_sample(sample, 10, 20)
    assert len(result) == 2
    np.testing.assert_array_equal(
        result["1"], np.array([3, 7, 5], dtype=np.int32)
    )
    np.testing.assert_array_equal(result["2"], np.array([30], dtype=np.int32))


def test_sex_correct_male():
    sample = {
        "1": np.array([1, 2]),
        "23": np.array([3, 4]),
        "24": np.array([5, 6]),
    }
    result = sex_correct(sample, Sex.MALE)
    np.testing.assert_array_equal(result["1"], np.array([1, 2]))
    np.testing.assert_array_equal(result["23"], np.array([6, 8]))
    np.testing.assert_array_equal(result["24"], np.array([10, 12]))


def test_sex_correct_female():
    sample = {
        "1": np.array([1, 2]),
        "23": np.array([3, 4]),
        "24": np.array([5, 6]),
    }
    result = sex_correct(sample, Sex.FEMALE)
    assert result is sample


def test_sex_correct_autosomal():
    sample = {
        "1": np.array([1, 2]),
        "23": np.array([3, 4]),
        "24": np.array([5, 6]),
    }
    result = sex_correct(sample, Sex.AUTOSOMAL)
    assert result is sample
