from __future__ import print_function, division

import struct

from astropy.tests.helper import pytest
import numpy as np

from numpy.testing import assert_array_almost_equal_nulp

from ..integrate import *


def almost_equal(a, b):
    c = struct.pack("<dd", a, b)
    d = struct.unpack("<qq", c)
    diff = abs(d[1] - d[0])
    return diff < 10

cases = [(1.0, 5.0), (1.0, 2.0), (1.0, 3.0), (1.5, 2.5), (1.2, 3.4),
         (1.6, 4.7), (2.2, 2.7), (1.0, 1.3), (1.0, 2.5), (1.3, 5.0),
         (3.3, 5.0), (4.4, 5.0), (1.0, 1.0), (5.0, 5.0), (2.0, 2.0),
         (2.0, 3.0), (2.0, 4.0), (2.0, 3.4), (3.0, 4.7), (3.0, 5.0)]


DTYPES = ['<f4', '>f4', '<f8', '>f8', '<i4', '>i4', '<i8', '>i8']


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_linear_subset(xmin, xmax):
    x = np.array([1., 2., 3., 4., 5.])
    y = 2. * x - 1.
    assert almost_equal(integrate_subset(x, y, xmin, xmax), (xmax * (xmax - 1) - xmin * (xmin - 1)))


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_loglog_subset(xmin, xmax):
    x = np.array([1., 2., 3., 4., 5.])
    y = 4. * x ** 3.
    assert almost_equal(integrate_loglog_subset(x, y, xmin, xmax), (xmax ** 4 - xmin ** 4))


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_loglog_subset_special1(xmin, xmax):
    # Special case for loglog is y = x^-1
    x = np.array([1., 2., 3., 4., 5.])
    y = 1. / x
    assert almost_equal(integrate_loglog_subset(x, y, xmin, xmax), np.log(xmax) - np.log(xmin))


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_loglog_subset_special2(xmin, xmax):
    # Special case for loglog is close to y = x^-1
    x = np.array([1., 2., 3., 4., 5.])
    y = x ** -1.1
    assert_array_almost_equal_nulp(integrate_loglog_subset(x, y, xmin, xmax), -10. * (xmax ** -0.1 - xmin ** -0.1), 50)


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_loglog_subset_special3(xmin, xmax):
    # Contains a zero value
    x = np.array([0., 2., 3., 4., 5.])
    y = 4. * x ** 3.
    assert almost_equal(integrate_loglog_subset(x, y, xmin, xmax), (max(xmax, 2.) ** 4 - max(xmin, 2.) ** 4))


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_loglin_subset(xmin, xmax):
    x = np.array([1., 2., 3., 4., 5.])
    y = np.log(2. * x)
    assert almost_equal(integrate_loglin_subset(x, y, xmin, xmax), (xmax * np.log(2. * xmax) - xmax - xmin * np.log(2. * xmin) + xmin))


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_linlog_subset(xmin, xmax):
    x = np.array([1., 2., 3., 4., 5.])
    y = np.exp(0.5 * x)
    assert almost_equal(integrate_linlog_subset(x, y, xmin, xmax), 2. * (np.exp(0.5 * xmax) - np.exp(0.5 * xmin)))


@pytest.mark.parametrize(('xmin', 'xmax'), cases)
def test_linlog_subset_special(xmin, xmax):
    # Special case for linlog is y[i] == y[i + 1]
    x = np.array([1., 2., 3., 4., 5.])
    y = np.repeat(1., x.shape)
    assert almost_equal(integrate_linlog_subset(x, y, xmin, xmax), xmax - xmin)


@pytest.mark.parametrize(('dtype_x', 'dtype_y'), zip(DTYPES, DTYPES))
def test_linear_types(dtype_x, dtype_y):
    x = np.array([1, 2], dtype=dtype_x)
    y = np.array([1, 1], dtype=dtype_y)
    assert integrate(x, y) == 1.


@pytest.mark.parametrize(('dtype_x', 'dtype_y'), zip(DTYPES, DTYPES))
def test_loglog_types(dtype_x, dtype_y):
    x = np.array([1, 2], dtype=dtype_x)
    y = np.array([1, 1], dtype=dtype_y)
    assert integrate_loglog(x, y) == 1.


@pytest.mark.parametrize(('dtype_x', 'dtype_y'), zip(DTYPES, DTYPES))
def test_loglin_types(dtype_x, dtype_y):
    x = np.array([1, 2], dtype=dtype_x)
    y = np.array([1, 1], dtype=dtype_y)
    assert integrate_loglin(x, y) == 1.


@pytest.mark.parametrize(('dtype_x', 'dtype_y'), zip(DTYPES, DTYPES))
def test_linlog_types(dtype_x, dtype_y):
    x = np.array([1, 2], dtype=dtype_x)
    y = np.array([1, 1], dtype=dtype_y)
    assert integrate_linlog(x, y) == 1.


def test_linear_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate(x, y)
    if isinstance(exc.value, basestring):
        assert exc.value == 'x is not monotonically increasing'
    else:
        assert exc.value.args[0] == 'x is not monotonically increasing'


def test_loglog_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate_loglog(x, y)
    if isinstance(exc.value, basestring):
        assert exc.value == 'x is not monotonically increasing'
    else:
        assert exc.value.args[0] == 'x is not monotonically increasing'


def test_loglin_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate_loglin(x, y)
    if isinstance(exc.value, basestring):
        assert exc.value == 'x is not monotonically increasing'
    else:
        assert exc.value.args[0] == 'x is not monotonically increasing'


def test_linlog_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate_linlog(x, y)
    if isinstance(exc.value, basestring):
        assert exc.value == 'x is not monotonically increasing'
    else:
        assert exc.value.args[0] == 'x is not monotonically increasing'
