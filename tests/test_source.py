import string
import random

import h5py
import atpy
import pytest
import numpy as np

from hyperion.sources import Source


def random_id(length=32):
    return string.join(random.sample(string.letters + string.digits, length), '')


def virtual_file():
    return h5py.File(random_id(), driver='core', backing_store=False)


def test_spectrum_atpy():
    t = atpy.Table()
    t.add_column('nu', [1, 2, 3])
    t.add_column('fnu', [1, 2, 3])
    s = Source()
    s.luminosity = 0.
    s.spectrum = t


def test_spectrum_atpy_invalid():
    t = atpy.Table()  # table is empty, so invalid
    s = Source()
    s.luminosity = 0.
    with pytest.raises(TypeError):
        s.spectrum = t


def test_spectrum_tuple():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    s.luminosity = 0.
    s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid1():
    nu = [1, 2, 3]  # should be a numpy array
    fnu = np.array([1, 2, 3])
    s = Source()
    s.luminosity = 0.
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid2():
    nu = np.array([1, 2, 3])
    fnu = [1, 2, 3]  # should be a numpy array
    s = Source()
    s.luminosity = 0.
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid3():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3, 4])  # sizes don't agree
    s = Source()
    s.luminosity = 0.
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid4():
    nu = np.array([[1, 2, 3], [4, 5, 6]])
    fnu = np.array([[1, 2, 3], [4, 5, 6]])  # arrays should be 1D
    s = Source()
    s.luminosity = 0.
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid5():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    s.luminosity = 0.
    with pytest.raises(TypeError):
        s.spectrum = (nu, fnu, fnu)  # too many items
