"""
Unit tests for the pure helpers in utils.py (no simulation required).
"""

import numpy as np
import pandas as pd
import pytest

from utils import stack, summarize, parse_unaids_value


def test_stack_assembles_expected_shape():
    years = np.arange(2000, 2005)
    dfs = [
        pd.DataFrame({'x': np.arange(len(years), dtype=float) + i})
        for i in range(3)
    ]
    out = stack(dfs, 'x', years)
    assert out.shape == (3, len(years))
    np.testing.assert_array_equal(out[0], np.arange(len(years)))


def test_stack_fills_nan_for_missing_column():
    years = np.arange(2000, 2005)
    dfs = [pd.DataFrame({'y': np.zeros(len(years))})]
    out = stack(dfs, 'missing', years)
    assert out.shape == (1, len(years))
    assert np.isnan(out).all()


def test_summarize_median_and_bands():
    arr = np.tile(np.arange(10, dtype=float), (5, 1))  # 5 identical seeds
    s = summarize(arr)
    np.testing.assert_array_equal(s['median'], np.arange(10))
    np.testing.assert_array_equal(s['lo'], np.arange(10))
    np.testing.assert_array_equal(s['hi'], np.arange(10))


def test_summarize_handles_nans():
    arr = np.array([[1.0, np.nan], [3.0, 5.0]])
    s = summarize(arr)
    assert s['median'][0] == 2.0
    assert s['median'][1] == 5.0  # only one non-nan value


@pytest.mark.parametrize('raw,expected', [
    ('1.2 [1.0 - 1.5]', 1.2),
    ('2900 [2400 - 3700]', 2900.0),
    ('<200', 200.0),
    ('11 000', 11000.0),
    ('...', None),
    ('', None),
])
def test_parse_unaids_value(raw, expected):
    result = parse_unaids_value(raw)
    if expected is None:
        assert np.isnan(result)
    else:
        assert result == expected


if __name__ == '__main__':
    pytest.main(['-x', '-v', __file__])
