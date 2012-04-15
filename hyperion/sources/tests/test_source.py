from __future__ import print_function, division

import string
import random

import h5py
import atpy
import pytest
import numpy as np

from .. import Source, PointSource, SpotSource, SphericalSource


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
    with pytest.raises(ValueError) as exc:
        s.write(v)
    assert exc.value.message == 'luminosity is not set'


def test_luminosity_invalid2():
    s = Source()
    with pytest.raises(ValueError) as exc:
        s.luminosity = np.array([1, 2, 3])  # luminosity should be a scalar
    assert exc.value.message == 'luminosity should be a scalar value'


def test_luminosity_invalid3():
    s = Source()
    with pytest.raises(ValueError) as exc:
        s.luminosity = 'invalid'  # luminosity should be a number
    assert exc.value.message == 'luminosity should be a numerical value'


def test_luminosity_invalid4():
    s = Source()
    with pytest.raises(ValueError) as exc:
        s.luminosity = -1.  # luminosity should be positive
    assert exc.value.message == 'luminosity should be positive'


def test_spectrum_atpy():
    t = atpy.Table()
    t.add_column('nu', [1, 2, 3])
    t.add_column('fnu', [1, 2, 3])
    s = Source()
    s.spectrum = t


def test_spectrum_atpy_invalid():
    t = atpy.Table()  # table is empty, so invalid
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.spectrum = t
    assert exc.value.message == 'spectrum ATpy Table does not contain a \'nu\' column'



def test_spectrum_tuple():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    s.spectrum = (nu, fnu)


def test_spectrum_tuple_invalid1():
    nu = [1, 2, 3]  # should be a numpy array
    fnu = np.array([1, 2, 3])
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.message == 'nu should be specified as a 1-D Numpy array'


def test_spectrum_tuple_invalid2():
    nu = np.array([1, 2, 3])
    fnu = [1, 2, 3]  # should be a numpy array
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.message == 'fnu should be specified as a 1-D Numpy array'



def test_spectrum_tuple_invalid3():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3, 4])  # sizes don't agree
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.message == 'nu and fnu should have the same shape'


def test_spectrum_tuple_invalid4():
    nu = np.array([[1, 2, 3], [4, 5, 6]])
    fnu = np.array([[1, 2, 3], [4, 5, 6]])  # arrays should be 1D
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.message == 'nu should be specified as a 1-D Numpy array'


def test_spectrum_tuple_invalid5():
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu, fnu)  # too many items
    assert exc.value.message == 'spectrum tuple or list should contain two elements'


# SpotSource

def test_point_longitude_none():
    s = SpotSource()
    s.longitude = None


def test_point_longitude_float():
    s = SpotSource()
    s.longitude = 1.


def test_spot_longitude_invalid1():
    v = virtual_file()
    s = SpotSource()
    # spot_longitude is not defined
    s.latitude = 1.
    s.radius = 1.
    s.temperature = 1.
    s.luminosity = 1.
    with pytest.raises(ValueError) as exc:
        s.write(v, 'test')
    assert exc.value.message == 'longitude is not set'


def test_spot_longitude_invalid2():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = np.array([1, 2, 3])  # longitude should be a scalar
    assert exc.value.message == 'longitude should be a scalar value'


def test_spot_longitude_invalid3():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = 'invalid'  # longitude should be a number
    assert exc.value.message == 'longitude should be a numerical value'


def test_spot_longitude_invalid4():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = -1.  # longitude should be in the range [0:360]
    assert exc.value.message == 'longitude should be in the range [0:360]'


def test_spot_longitude_invalid5():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = 366.  # longitude should be in the range [0:360]
    assert exc.value.message == 'longitude should be in the range [0:360]'


def test_spot_latitude_none():
    s = SpotSource()
    s.latitude = None


def test_spot_latitude_float1():
    s = SpotSource()
    s.latitude = -80.


def test_spot_latitude_float2():
    s = SpotSource()
    s.latitude = 70.


def test_spot_latitude_invalid1():
    v = virtual_file()
    s = SpotSource()
    s.longitude = 1.
    # spot_latitude is not defined
    s.radius = 1.
    s.temperature = 1.
    s.luminosity = 1.
    with pytest.raises(ValueError) as exc:
        s.write(v, 'test')
    assert exc.value.message == 'latitude is not set'


def test_spot_latitude_invalid2():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = np.array([1, 2, 3])  # latitude should be a scalar
    assert exc.value.message == 'latitude should be a scalar value'


def test_spot_latitude_invalid3():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = 'invalid'  # latitude should be a number
    assert exc.value.message == 'latitude should be a numerical value'


def test_spot_latitude_invalid4():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = -92.  # latitude should be in the range [-90:90]
    assert exc.value.message == 'latitude should be in the range [-90:90]'


def test_spot_latitude_invalid5():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = 100  # latitude should be in the range [-90:90]
    assert exc.value.message == 'latitude should be in the range [-90:90]'


def test_spot_radius_none():
    s = SpotSource()
    s.radius = None


def test_spot_radius_float():
    s = SpotSource()
    s.radius = 1.e10


def test_spot_radius_invalid1():
    v = virtual_file()
    s = SpotSource()
    s.longitude = 1.
    s.latitude = 1.
    # spot_radius is not defined
    s.temperature = 1.
    s.luminosity = 1.
    with pytest.raises(ValueError) as exc:
        s.write(v, 'test')
    assert exc.value.message == 'radius is not set'


def test_spot_radius_invalid2():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.radius = np.array([1, 2, 3])  # radius should be a scalar
    assert exc.value.message == 'radius should be a scalar value'


def test_spot_radius_invalid3():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.radius = 'invalid'  # radius should be a number
    assert exc.value.message == 'radius should be a numerical value'


def test_spot_radius_invalid4():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.radius = -1.  # radius should be positive
    assert exc.value.message == 'radius should be positive'


# PointSource

def test_point_spectrum_tuple_invalid():
    # This is to check that the __setattr__ for Source still gets called after
    # inheritance.
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu, fnu)  # too many items
    assert exc.value.message == 'spectrum tuple or list should contain two elements'


def test_point_position_none():
    s = PointSource()
    s.position = None


def test_point_position_tuple():
    s = PointSource()
    s.position = (0., 1., 2.)


def test_point_position_tuple_invalid():
    s = PointSource()
    with pytest.raises(ValueError) as exc:
        s.position = (0., 1., 2., 4.)  # too many elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_point_position_list():
    s = PointSource()
    s.position = [1., 2., 3.]


def test_point_position_list_invalid():
    s = PointSource()
    with pytest.raises(ValueError) as exc:
        s.position = [1., 2.]  # too few elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_point_position_numpy():
    s = PointSource()
    s.position = np.array([2., 3., 4.])


def test_point_position_numpy_invalid1():
    s = PointSource()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([2.])  # too few elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_point_position_numpy_invalid2():
    s = PointSource()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([[1., 2., 3.]])  # wrong dimensionality
    assert exc.value.message == 'position should be a 1-D sequence'

# SphericalSource


def test_spherical_position_none():
    s = SphericalSource()
    s.position = None


def test_spherical_position_tuple():
    s = SphericalSource()
    s.position = (0., 1., 2.)


def test_spherical_position_tuple_invalid():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = (0., 1., 2., 4.)  # too many elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_spherical_position_list():
    s = SphericalSource()
    s.position = [1., 2., 3.]


def test_spherical_position_list_invalid():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = [1., 2.]  # too few elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_spherical_position_numpy():
    s = SphericalSource()
    s.position = np.array([2., 3., 4.])


def test_spherical_position_numpy_invalid1():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([2.])  # too few elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_spherical_position_numpy_invalid2():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([[1., 2., 3.]])  # wrong dimensionality
    assert exc.value.message == 'position should be a 1-D sequence'


def test_spherical_radius_none():
    s = SphericalSource()
    s.radius = None


def test_spherical_radius_float():
    s = SphericalSource()
    s.radius = 1.e10


def test_spherical_radius_invalid1():
    v = virtual_file()
    s = SphericalSource()
    s.position = (0., 0., 0.)
    # radius is not defined
    s.limb = True
    s.temperature = 1.
    s.luminosity = 1.
    with pytest.raises(ValueError) as exc:
        s.write(v, 'test')
    assert exc.value.message == 'radius is not set'


def test_spherical_radius_invalid2():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.radius = np.array([1, 2, 3])  # radius should be a scalar
    assert exc.value.message == 'radius should be a scalar value'


def test_spherical_radius_invalid3():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.radius = 'invalid'  # radius should be a number
    assert exc.value.message == 'radius should be a numerical value'


def test_spherical_radius_invalid4():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.radius = -1.  # radius should be positive
    assert exc.value.message == 'radius should be positive'


def test_spherical_limb_none():
    s = SphericalSource()
    s.limb = None


def test_spherical_limb_true():
    s = SphericalSource()
    s.limb = True


def test_spherical_limb_false():
    s = SphericalSource()
    s.limb = False


def test_spherical_limb_invalid2():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.limb = np.array([1, 2, 3])  # limb should be a boolean
    assert exc.value.message == 'limb should be a boolean value (True/False)'


def test_spherical_limb_invalid3():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.limb = 'invalid'  # limb should be a boolean
    assert exc.value.message == 'limb should be a boolean value (True/False)'


def test_spherical_limb_invalid4():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.limb = 1  # limb should be a boolean
    assert exc.value.message == 'limb should be a boolean value (True/False)'


def test_spherical_limb_invalid5():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.limb = 1.  # limb should be a boolean
    assert exc.value.message == 'limb should be a boolean value (True/False)'


# ExternalSpherical

def test_external_spherical_position_none():
    s = SphericalSource()
    s.position = None


def test_external_spherical_position_tuple():
    s = SphericalSource()
    s.position = (0., 1., 2.)


def test_external_spherical_position_tuple_invalid():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = (0., 1., 2., 4.)  # too many elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_external_spherical_position_list():
    s = SphericalSource()
    s.position = [1., 2., 3.]


def test_external_spherical_position_list_invalid():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = [1., 2.]  # too few elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_external_spherical_position_numpy():
    s = SphericalSource()
    s.position = np.array([2., 3., 4.])


def test_external_spherical_position_numpy_invalid1():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([2.])  # too few elements
    assert exc.value.message == 'position should be a sequence of 3 values'


def test_external_spherical_position_numpy_invalid2():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([[1., 2., 3.]])  # wrong dimensionality
    assert exc.value.message == 'position should be a 1-D sequence'


def test_external_spherical_radius_none():
    s = SphericalSource()
    s.radius = None


def test_external_spherical_radius_float():
    s = SphericalSource()
    s.radius = 1.e10


def test_external_spherical_radius_invalid1():
    v = virtual_file()
    s = SphericalSource()
    s.position = (0., 0., 0.)
    # radius is not defined
    s.limb = True
    s.temperature = 1.
    s.luminosity = 1.
    with pytest.raises(ValueError) as exc:
        s.write(v, 'test')
    assert exc.value.message == 'radius is not set'


def test_external_spherical_radius_invalid2():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.radius = np.array([1, 2, 3])  # radius should be a scalar
    assert exc.value.message == 'radius should be a scalar value'


def test_external_spherical_radius_invalid3():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.radius = 'invalid'  # radius should be a number
    assert exc.value.message == 'radius should be a numerical value'


def test_external_spherical_radius_invalid4():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.radius = -1.  # radius should be positive
    assert exc.value.message == 'radius should be positive'
