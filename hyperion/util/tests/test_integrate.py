from __future__ import print_function, division

import struct

import pytest
import numpy as np

from ..integrate import *


def almost_equal(a, b):
    c = struct.pack("<dd", a, b)
    d = struct.unpack("<qq", c)
    diff = abs(d[1] - d[0])
    return diff < 10

cases = [(1.0, 5.0), (1.0, 2.0), (1.0, 3.0), (1.5, 2.5), (1.2, 3.4),
         (1.6, 4.7), (2.2, 2.7), (1.0, 1.3), (1.0, 2.5), (1.3, 5.0),
         (3.3, 5.0), (4.4, 5.0), (1.0, 1.0), (5.0, 5.0), (2.0, 2.0)]


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
def test_loglog_subset_special(xmin, xmax):
    # Special case for loglog is y = x^-1
    x = np.array([1., 2., 3., 4., 5.])
    y = 1. / x
    assert almost_equal(integrate_loglog_subset(x, y, xmin, xmax), np.log(xmax) - np.log(xmin))


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


def test_linear_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate(x, y)
    assert exc.value.args[0] == 'x is not monotonically increasing'


def test_loglog_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate_loglog(x, y)
    assert exc.value.args[0] == 'x is not monotonically increasing'


def test_loglin_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate_loglin(x, y)
    assert exc.value.args[0] == 'x is not monotonically increasing'


def test_linlog_not_monotonic():
    x = np.array([1., 2., 4., 3., 5.])
    y = 2. * x - 1.
    with pytest.raises(ValueError) as exc:
        integrate_linlog(x, y)
    assert exc.value.args[0] == 'x is not monotonically increasing'
