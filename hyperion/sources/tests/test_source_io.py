from __future__ import print_function, division

import numpy as np
from numpy.testing import assert_equal
from astropy.tests.helper import pytest

from .. import (Source, PointSource, PointSourceCollection, SpotSource,
                SphericalSource, ExternalSphericalSource, ExternalBoxSource,
                MapSource, PlaneParallelSource, read_source)
from ...grid import CartesianGrid
from ...util.functions import virtual_file


def test_io_source_temperature():
    s1 = Source()
    s1.luminosity = 1.
    s1.temperature = 5000.
    v = virtual_file()
    s1.write(v)
    s2 = Source.read(v)
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None


def test_io_source_spectrum():
    s1 = Source()
    s1.luminosity = 1.
    s1.spectrum = ([1., 2., 3.], [4., 5., 6.])
    v = virtual_file()
    s1.write(v)
    s2 = Source.read(v)
    assert s2.luminosity == s1.luminosity
    assert s2.temperature is None
    assert_equal(s2.spectrum['nu'], s1.spectrum['nu'])
    assert_equal(s2.spectrum['fnu'], s1.spectrum['fnu'])


def test_io_source_lte():
    s1 = Source()
    s1.luminosity = 1.
    v = virtual_file()
    s1.write(v)
    s2 = Source.read(v)
    assert s2.luminosity == s1.luminosity
    assert s2.temperature is None
    assert s2.spectrum is None


def test_io_spot_source():
    s1 = SpotSource()
    s1.luminosity = 1.
    s1.temperature = 5000.
    s1.longitude = 12.
    s1.latitude = 34.
    s1.radius = 5.
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert s2.longitude == s1.longitude
    assert s2.latitude == s1.latitude
    assert s2.radius == s1.radius


def test_io_point_source():
    s1 = PointSource()
    s1.luminosity = 1.
    s1.temperature = 5000.
    s1.position = [1., 2., 3.]
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert_equal(s2.position, s1.position)


def test_io_point_source_collection():
    s1 = PointSourceCollection()
    s1.luminosity = np.array([1.,3.,4.])
    s1.temperature = 5000.
    s1.position = np.array([[3., 2., 5.], [-3., 2., 6.], [9., 2., 1.]])
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert_equal(s2.luminosity, s1.luminosity)
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert_equal(s2.position, s1.position)

@pytest.mark.parametrize('limb', [True, False])
def test_io_spherical_source(limb):
    s1 = SphericalSource()
    s1.luminosity = 1.
    s1.temperature = 5000.
    s1.position = [1., 2., 3.]
    s1.radius = 4.
    s1.limb = limb
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert_equal(s2.position, s1.position)
    assert s2.radius == s1.radius
    assert s2.limb is s1.limb


def test_io_spherical_source_with_spots():
    s1 = SphericalSource()
    s1.luminosity = 1.
    s1.temperature = 5000.
    s1.position = [1., 2., 3.]
    s1.radius = 4.
    s1.limb = True
    spot1 = s1.add_spot()
    spot1.luminosity = 0.5
    spot1.temperature = 3000.
    spot1.longitude = 1.
    spot1.latitude = 2.
    spot1.radius = 3.
    spot2 = s1.add_spot()
    spot2.luminosity = 0.6
    spot2.temperature = 2000.
    spot2.longitude = 4.
    spot2.latitude = 5.
    spot2.radius = 6.
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert_equal(s2.position, s1.position)
    assert s2.radius == s1.radius
    assert s2.limb is s1.limb
    for ispot in range(2):
        assert s2._spots[ispot].longitude == s1._spots[ispot].longitude
        assert s2._spots[ispot].latitude == s1._spots[ispot].latitude
        assert s2._spots[ispot].radius == s1._spots[ispot].radius


def test_io_external_spherical_source():
    s1 = ExternalSphericalSource()
    s1.luminosity = 1.
    s1.temperature = 5000.
    s1.position = [1., 2., 3.]
    s1.radius = 4.
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert_equal(s2.position, s1.position)
    assert s2.radius == s1.radius


def test_io_external_box_source():
    s1 = ExternalBoxSource()
    s1.luminosity = 1.
    s1.temperature = 5000.
    s1.bounds = [(1., 2.), (3., 4.), (5., 6.)]
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert_equal(s2.bounds, s1.bounds)


def test_io_map_source():

    # Set up coordinate grid
    g = CartesianGrid([-1., 0., 1.], [-1., 1.], [-1., -0.2, 0.2, 1.])

    s1 = MapSource()
    s1.luminosity = 1.
    s1.map = np.ones((3, 1, 2))
    v = virtual_file()
    s1.write(v, 'test', g)

    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert s2.luminosity == s1.luminosity
    assert_equal(s2.map, s1.map)


def test_io_plane_parallel_source():
    s1 = PlaneParallelSource()
    s1.luminosity = 1.
    s1.temperature = 5000.
    s1.position = [1., 2., 3.]
    s1.radius = 4.
    s1.direction = (5., 6.)
    v = virtual_file()
    s1.write(v, 'test')
    s2 = read_source(v['test'])
    assert s2.name == s1.name
    assert s2.luminosity == s1.luminosity
    assert s2.temperature == s1.temperature
    assert s2.spectrum is None
    assert_equal(s2.position, s1.position)
    assert s2.radius == s1.radius
    assert_equal(s2.direction, s1.direction)
