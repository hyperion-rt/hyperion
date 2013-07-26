from astropy.tests.helper import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp

from ..interpolate import *


DTYPES = ['<f4', '>f4', '<f8', '>f8', '<i4', '>i4', '<i8', '>i8']


class GenericTests(object):

    def setup_class(self):
        raise Exception("This class should not be used directly")

    def test_interp_linear_array(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = np.linspace(1., 10., 100)
        assert_array_almost_equal_nulp(self.interp(x, y, xval), self.f(xval), 15)

    def test_interp_linear_scalar1(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = 1.
        assert_array_almost_equal_nulp(self.interp(x, y, xval), self.f(xval), 15)

    def test_interp_linear_scalar2(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = 4.
        assert_array_almost_equal_nulp(self.interp(x, y, xval), self.f(xval), 15)

    def test_interp_linear_scalar3(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = 10.
        assert_array_almost_equal_nulp(self.interp(x, y, xval), self.f(xval), 15)

    def test_interp_linear_array_invalid1(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = np.linspace(0., 7., 100)
        with pytest.raises(Exception) as exc:
            self.interp(x, y, xval)
        assert exc.value.args[0] == 'x values are out of interpolation bounds'

    def test_interp_linear_array_invalid2(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = np.linspace(2., 11., 100)
        with pytest.raises(Exception) as exc:
            self.interp(x, y, xval)
        assert exc.value.args[0] == 'x values are out of interpolation bounds'

    def test_interp_linear_scalar_invalid1(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = 0.9
        with pytest.raises(Exception) as exc:
            self.interp(x, y, xval)
        assert exc.value.args[0] == 'x value is out of interpolation bounds'

    def test_interp_linear_scalar_invalid2(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = 11.
        with pytest.raises(Exception) as exc:
            self.interp(x, y, xval)
        assert exc.value.args[0] == 'x value is out of interpolation bounds'

    def test_interp_linear_array_fill(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = np.linspace(0., 10., 100)
        ref = self.f(xval)
        ref[xval < x[0]] = -2.
        assert_array_almost_equal_nulp(self.interp(x, y, xval, bounds_error=False, fill_value=-2.), ref, 15)

    def test_interp_linear_scalar_fill(self):
        x = np.linspace(1., 10., 10)
        y = self.f(x)
        xval = 11.
        ref = -2.
        assert_array_almost_equal_nulp(self.interp(x, y, xval, bounds_error=False, fill_value=-2.), ref, 15)

    def test_length_mismatch(self):
        x = np.linspace(1., 10., 10)
        y = np.linspace(1., 10., 11)
        xval = np.linspace(1., 10., 100)
        with pytest.raises(Exception) as exc:
            self.interp(x, y, xval)
        assert exc.value.args[0] == 'x and y should have the same length'

    @pytest.mark.parametrize(('dtype_x', 'dtype_y'), zip(DTYPES, DTYPES))
    def test_types(self, dtype_x, dtype_y):
        x = np.array([1, 5], dtype=dtype_x)
        y = np.array([1, 1], dtype=dtype_y)
        xval = 3.
        assert self.interp(x, y, xval) == 1.


class TestLinear(GenericTests):

    def setup_class(self):
        self.f = lambda self, x: 2 * x - 1
        self.interp = interp1d_fast


class TestLogLog(GenericTests):

    def setup_class(self):
        self.f = lambda self, x: 4. * x ** 3
        self.interp = interp1d_fast_loglog


class TestLogLin(GenericTests):

    def setup_class(self):
        self.f = lambda self, x: np.log(2. * x)
        self.interp = interp1d_fast_loglin


class TestLinLog(GenericTests):

    def setup_class(self):
        self.f = lambda self, x: np.exp(0.5 * x)
        self.interp = interp1d_fast_linlog
