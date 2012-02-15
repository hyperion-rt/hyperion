import string
import random

import h5py
import atpy
import pytest
import numpy as np

from hyperion.sources import Source, PointSource


def random_id(length=32):
    return string.join(random.sample(string.letters + string.digits, length), '')


def virtual_file():
    return h5py.File(random_id(), driver='core', backing_store=False)


# Source


def test_luminosity():
    v = virtual_file()
    s = Source()
    s.luminosity = 1.
    s.write(v)


def test_luminosity_invalid1():
    v = virtual_file()
    s = Source()
    # luminosity is not defined
    with pytest.raises(ValueError):
        s.write(v)


def test_luminosity_invalid2():
    v = virtual_file()
    s = Source()
    s.luminosity = np.array([1, 2, 3])  # luminosity should be a scalar
    with pytest.raises(ValueError):
        s.write(v)


def test_luminosity_invalid3():
    v = virtual_file()
    s = Source()
    s.luminosity = 'invalid'  # luminosity should be a number
    with pytest.raises(ValueError):
        s.write(v)


def test_spectrum_atpy():
    t = atpy.Table()
    t.add_column('nu', [1, 2, 3])
    t.add_column('fnu', [1, 2, 3])
    s = Source()
    s.spectrum = t


def test_spectrum_atpy_invalid():
    t = atpy.Table()  # table is empty, so invalid
    s = Source()
    with pytest.raises(TypeError):
        s.spectrum = t


def test_spectrum_tuple():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid1():
    nu = [1, 2, 3]  # should be a numpy array
    fnu = np.array([1, 2, 3])
    s = Source()
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid2():
    nu = np.array([1, 2, 3])
    fnu = [1, 2, 3]  # should be a numpy array
    s = Source()
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid3():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3, 4])  # sizes don't agree
    s = Source()
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid4():
    nu = np.array([[1, 2, 3], [4, 5, 6]])
    fnu = np.array([[1, 2, 3], [4, 5, 6]])  # arrays should be 1D
    s = Source()
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid5():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu, fnu)  # too many items


# PointSource

def test_point_spectrum_tuple_invalid():
    # This is to check that the __setattr__ for Source still gets called after
    # inheritance.
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu, fnu)  # too many items


def test_point_position_none():
    s = PointSource()
    s.position = None


def test_point_position_tuple():
    s = PointSource()
    s.position = (0., 1., 2.)


def test_point_position_tuple_invalid():
    s = PointSource()
    with pytest.raises(ValueError):
        s.position = (0., 1., 2., 4.)  # too many elements


def test_point_position_list():
    s = PointSource()
    s.position = [1., 2., 3.]


def test_point_position_list_invalid():
    s = PointSource()
    with pytest.raises(ValueError):
        s.position = [1., 2.]  # too few elements


def test_point_position_numpy():
    s = PointSource()
    s.position = np.array([2., 3., 4.])


def test_point_position_numpy_invalid1():
    s = PointSource()
    with pytest.raises(ValueError):
        s.position = np.array([2.])  # too few elements


def test_point_position_numpy_invalid2():
    s = PointSource()
    with pytest.raises(ValueError):
        s.position = np.array([[1., 2., 3.]])  # wrong dimensionality
