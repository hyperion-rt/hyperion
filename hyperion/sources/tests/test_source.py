from __future__ import print_function, division

from astropy.tests.helper import pytest
import numpy as np
from astropy.table import Table, Column

from .. import (Source, PointSource, PointSourceCollection, SpotSource,
                SphericalSource, ExternalSphericalSource, ExternalBoxSource,
                MapSource, PlaneParallelSource)

from ...grid import (CartesianGrid,
                    CylindricalPolarGrid,
                    SphericalPolarGrid,
                    AMRGrid,
                    OctreeGrid)

from ...util.functions import virtual_file

ALL_SOURCES = [Source, PointSource, PointSourceCollection, SpotSource,
               SphericalSource, ExternalSphericalSource, ExternalBoxSource,
               MapSource, PlaneParallelSource]


# SCALAR LUMINOSITY


@pytest.mark.parametrize(('source_type'), list(set(ALL_SOURCES) - set([PointSourceCollection])))
def test_luminosity_scalar(source_type):
    s = source_type()
    s.luminosity = 1.

@pytest.mark.parametrize(('source_type'), list(set(ALL_SOURCES) - set([PointSourceCollection])))
def test_luminosity_scalar_invalid2(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.luminosity = np.array([1, 2, 3])  # luminosity should be a scalar
    assert exc.value.args[0] == 'luminosity should be a scalar value'


@pytest.mark.parametrize(('source_type'), list(set(ALL_SOURCES) - set([PointSourceCollection])))
def test_luminosity_scalar_invalid3(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.luminosity = 'invalid'  # luminosity should be a number
    assert exc.value.args[0] == 'luminosity should be a numerical value'


@pytest.mark.parametrize(('source_type'), list(set(ALL_SOURCES) - set([PointSourceCollection])))
def test_luminosity_scalar_invalid4(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.luminosity = -1.  # luminosity should be positive
    assert exc.value.args[0] == 'luminosity should be positive'

# ARRAY LUMINOSITY

@pytest.mark.parametrize(('source_type'), [PointSourceCollection])
def test_luminosity_array(source_type):
    s = source_type()
    s.luminosity = np.array([1, 2, 3])

@pytest.mark.parametrize(('source_type'), [PointSourceCollection])
def test_luminosity_array_invalid1(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.luminosity = 1.
    assert exc.value.args[0] == 'luminosity should be a Numpy array'

@pytest.mark.parametrize(('source_type'), [PointSourceCollection])
def test_luminosity_array_invalid2(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.luminosity = np.ones((2,2))
    assert exc.value.args[0] == 'luminosity should be a 1-D array'

@pytest.mark.parametrize(('source_type'), [PointSourceCollection])
def test_luminosity_array_invalid3(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.luminosity = np.array([1, 2, -3])
    assert exc.value.args[0] == 'luminosity should be positive'

# TEMPERATURE


@pytest.mark.parametrize(('source_type'), list(set(ALL_SOURCES) - set([MapSource])))
def test_temperature(source_type):
    v = virtual_file()
    s = source_type()
    s.temperature = 1.


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_temperature_invalid2(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.temperature = np.array([1, 2, 3])  # temperature should be a scalar
    assert exc.value.args[0] == 'temperature should be a scalar value'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_temperature_invalid3(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.temperature = 'invalid'  # temperature should be a number
    assert exc.value.args[0] == 'temperature should be a numerical value'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_temperature_invalid4(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.temperature = -1.  # temperature should be positive
    assert exc.value.args[0] == 'temperature should be positive'

# SPECTRUM


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_astropy(source_type):
    t = Table()
    t.add_column(Column(data=[1, 2, 3], name='nu'))
    t.add_column(Column(data=[1, 2, 3], name='fnu'))
    s = source_type()
    s.spectrum = t


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_astropy_invalid(source_type):
    t = Table()  # table is empty, so invalid
    s = source_type()
    with pytest.raises(TypeError) as exc:
        s.spectrum = t
    assert exc.value.args[0] == 'spectrum Table does not contain a \'nu\' column'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_valid1(source_type):
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = source_type()
    s.spectrum = (nu, fnu)


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_valid2(source_type):
    nu = [1, 2, 3]
    fnu = np.array([1, 2, 3])
    s = source_type()
    s.spectrum = (nu, fnu)


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_valid3(source_type):
    nu = np.array([1, 2, 3])
    fnu = [1, 2, 3]
    s = source_type()
    s.spectrum = (nu, fnu)


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_sort(source_type):
    nu = [3, 1, 2]
    fnu = [1, 2, 3]
    s = source_type()
    s.spectrum = (nu, fnu)
    assert np.all(s.spectrum['nu'] == np.array([1, 2, 3]))
    assert np.all(s.spectrum['fnu'] == np.array([2, 3, 1]))


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_invalid1(source_type):
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3, 4])  # sizes don't agree
    s = source_type()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.args[0] == 'nu and fnu should have the same shape'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_invalid2(source_type):
    nu = np.array([[1, 2, 3], [4, 5, 6]])
    fnu = np.array([[1, 2, 3], [4, 5, 6]])  # arrays should be 1D
    s = source_type()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.args[0] == 'nu should be a 1-D sequence'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_invalid3(source_type):
    nu = np.array([1, 2, 3])
    fnu = np.array([1, 2, 3])
    s = source_type()
    with pytest.raises(TypeError) as exc:
        s.spectrum = (nu, fnu, fnu)  # too many items
    assert exc.value.args[
        0] == 'spectrum tuple or list should contain two elements'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_tuple_invalid4(source_type):
    nu = np.array([1, 2, 2])  # duplicate values
    fnu = np.array([1, 2, 3])
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.args[0] == 'nu sequence contains duplicate values'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_set_temperature_spectrum(source_type):
    s = source_type()
    s.temperature = 1000.
    with pytest.raises(Exception) as exc:
        s.spectrum = (np.array([1, 2, 3]), np.array([4, 5, 6]))  # temperature has already been specified
    assert exc.value.args[0] == 'A temperature has already been set, so cannot set a spectrum'


@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_set_spectrum_temperature(source_type):
    s = source_type()
    s.spectrum = (np.array([1, 2, 3]), np.array([4, 5, 6]))
    with pytest.raises(Exception) as exc:
        s.temperature = 1000.  # spectrum has already been specified
    assert exc.value.args[0] == 'A spectrum has already been set, so cannot set a temperature'

@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_negative(source_type):

    nu = np.array([1, 2, -3])
    fnu = np.array([1, 2, 3])
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.args[0] == 'nu should be strictly positive'

    nu = np.array([1, 2, 3])
    fnu = np.array([1, -2, 3])
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.args[0] == 'fnu should be positive'

@pytest.mark.parametrize(('source_type'), ALL_SOURCES)
def test_spectrum_nan(source_type):

    nu = np.array([1, 2, np.nan])
    fnu = np.array([1, 2, 3])
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.args[0] == 'nu contains NaN/Inf values'

    nu = np.array([1, 2, 3])
    fnu = np.array([1, np.inf, 3])
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.spectrum = (nu, fnu)
    assert exc.value.args[0] == 'fnu contains NaN/Inf values'

# POSITION

SOURCES_POSITION = [PointSource,
                    SphericalSource,
                    ExternalSphericalSource,
                    PlaneParallelSource]


def test_position_tests_complete():
    expected = [x for x in ALL_SOURCES if hasattr(x, 'position')]
    extra = list(set(SOURCES_POSITION) - set(expected))
    assert extra == []


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_none(source_type):
    s = source_type()
    s.position = None


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_tuple(source_type):
    s = source_type()
    s.position = (0., 1., 2.)


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_tuple_invalid(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.position = (0., 1., 2., 4.)  # too many elements
    assert exc.value.args[0] == 'position should be a sequence of 3 values'


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_list(source_type):
    s = source_type()
    s.position = [1., 2., 3.]


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_list_invalid(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.position = [1., 2.]  # too few elements
    assert exc.value.args[0] == 'position should be a sequence of 3 values'


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_numpy(source_type):
    s = source_type()
    s.position = np.array([2., 3., 4.])


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_numpy_invalid1(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([2.])  # too few elements
    assert exc.value.args[0] == 'position should be a sequence of 3 values'


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION)
def test_position_numpy_invalid2(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([[1., 2., 3.]])  # wrong dimensionality
    assert exc.value.args[0] == 'position should be a 1-D sequence'

# POSITION (ARRAY)

SOURCES_POSITION_ARRAY = [PointSourceCollection]


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION_ARRAY)
def test_position_none(source_type):
    s = source_type()
    s.position = None


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION_ARRAY)
def test_position_array(source_type):
    s = source_type()
    s.position = np.array([[1., 2., 3.]])


@pytest.mark.parametrize(('source_type'), SOURCES_POSITION_ARRAY)
def test_position_invalid1(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.position = [0., 1., 2.]
    assert exc.value.args[0] == 'position should be a Numpy array'

@pytest.mark.parametrize(('source_type'), SOURCES_POSITION_ARRAY)
def test_position_invalid2(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([0., 1., 2.])
    assert exc.value.args[0] == 'position should be a 2-D array'

@pytest.mark.parametrize(('source_type'), SOURCES_POSITION_ARRAY)
def test_position_invalid3(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([[0., 1., 2., 4.]])
    assert exc.value.args[0] == 'position should be an Nx3 array'

# RADIUS

SOURCES_RADIUS = [SphericalSource,
                  ExternalSphericalSource,
                  PlaneParallelSource]


def test_radius_tests_complete():
    expected = [x for x in ALL_SOURCES if hasattr(x, 'radius')]
    extra = list(set(SOURCES_RADIUS) - set(expected))
    assert extra == []


@pytest.mark.parametrize(('source_type'), SOURCES_RADIUS)
def test_radius_none(source_type):
    s = source_type()
    s.radius = None


@pytest.mark.parametrize(('source_type'), SOURCES_RADIUS)
def test_radius_float(source_type):
    s = source_type()
    s.radius = 1.e10


@pytest.mark.parametrize(('source_type'), SOURCES_RADIUS)
def test_radius_invalid2(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.radius = np.array([1, 2, 3])  # radius should be a scalar
    assert exc.value.args[0] == 'radius should be a scalar value'


@pytest.mark.parametrize(('source_type'), SOURCES_RADIUS)
def test_radius_invalid3(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.radius = 'invalid'  # radius should be a number
    assert exc.value.args[0] == 'radius should be a numerical value'


@pytest.mark.parametrize(('source_type'), SOURCES_RADIUS)
def test_radius_invalid4(source_type):
    s = source_type()
    with pytest.raises(ValueError) as exc:
        s.radius = -1.  # radius should be positive
    assert exc.value.args[0] == 'radius should be positive'

# PointSource

REQUIRED = {}

REQUIRED[Source] = {}
REQUIRED[Source]['luminosity'] = 1.

REQUIRED[PointSource] = {}
REQUIRED[PointSource]['luminosity'] = 1.
REQUIRED[PointSource]['temperature'] = 1.

REQUIRED[SphericalSource] = {}
REQUIRED[SphericalSource]['luminosity'] = 1.
REQUIRED[SphericalSource]['temperature'] = 1.
REQUIRED[SphericalSource]['radius'] = 1.

REQUIRED[ExternalSphericalSource] = {}
REQUIRED[ExternalSphericalSource]['luminosity'] = 1.
REQUIRED[ExternalSphericalSource]['temperature'] = 1.
REQUIRED[ExternalSphericalSource]['radius'] = 1.

REQUIRED[ExternalBoxSource] = {}
REQUIRED[ExternalBoxSource]['luminosity'] = 1.
REQUIRED[ExternalBoxSource]['temperature'] = 1.
REQUIRED[ExternalBoxSource]['bounds'] = [[1., 2.], [3., 4.], [5., 6.]]

REQUIRED[PlaneParallelSource] = {}
REQUIRED[PlaneParallelSource]['luminosity'] = 1.
REQUIRED[PlaneParallelSource]['temperature'] = 1.
REQUIRED[PlaneParallelSource]['radius'] = 1.
REQUIRED[PlaneParallelSource]['direction'] = (1., 2.)

# Note that position is not required since it defaults to the origin.
# limb is also not required, since it defaults to False.

# Test that no errors are raised if all attributes are present


@pytest.mark.parametrize(('source_type'), REQUIRED)
def test_all(source_type):
    v = virtual_file()
    s = source_type()
    for attribute in REQUIRED[source_type]:
        setattr(s, attribute, REQUIRED[source_type][attribute])
    if source_type == Source:
        s.write(v)
    else:
        s.write(v, 'test')

# Test that an error is raised if one attribute is missing

combinations = [(s, a) for s in REQUIRED for a in REQUIRED[s]]


@pytest.mark.parametrize(('source_type', 'missing'), combinations)
def test_missing(source_type, missing):
    v = virtual_file()
    s = source_type()
    for attribute in REQUIRED[source_type]:
        if attribute == missing:
            continue
        setattr(s, attribute, REQUIRED[source_type][attribute])
    with pytest.raises(ValueError) as exc:
        if source_type == Source:
            s.write(v)
        else:
            s.write(v, 'test')
    if missing == 'temperature':
        assert exc.value.args[0].endswith('cannot have LTE spectrum')
    else:
        assert exc.value.args[0] == '{0:s} is not set'.format(missing)

# SpotSource - specific tests


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
    assert exc.value.args[0] == 'longitude is not set'


def test_spot_longitude_invalid2():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = np.array([1, 2, 3])  # longitude should be a scalar
    assert exc.value.args[0] == 'longitude should be a scalar value'


def test_spot_longitude_invalid3():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = 'invalid'  # longitude should be a number
    assert exc.value.args[0] == 'longitude should be a numerical value'


def test_spot_longitude_invalid4():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = -1.  # longitude should be in the range [0:360]
    assert exc.value.args[0] == 'longitude should be in the range [0:360]'


def test_spot_longitude_invalid5():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.longitude = 366.  # longitude should be in the range [0:360]
    assert exc.value.args[0] == 'longitude should be in the range [0:360]'


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
    assert exc.value.args[0] == 'latitude is not set'


def test_spot_latitude_invalid2():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = np.array([1, 2, 3])  # latitude should be a scalar
    assert exc.value.args[0] == 'latitude should be a scalar value'


def test_spot_latitude_invalid3():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = 'invalid'  # latitude should be a number
    assert exc.value.args[0] == 'latitude should be a numerical value'


def test_spot_latitude_invalid4():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = -92.  # latitude should be in the range [-90:90]
    assert exc.value.args[0] == 'latitude should be in the range [-90:90]'


def test_spot_latitude_invalid5():
    s = SpotSource()
    with pytest.raises(ValueError) as exc:
        s.latitude = 100  # latitude should be in the range [-90:90]
    assert exc.value.args[0] == 'latitude should be in the range [-90:90]'

# SphericalSource - specific tests


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
    assert exc.value.args[0] == 'limb should be a boolean value (True/False)'


def test_spherical_limb_invalid3():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.limb = 'invalid'  # limb should be a boolean
    assert exc.value.args[0] == 'limb should be a boolean value (True/False)'


def test_spherical_limb_invalid4():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.limb = 1  # limb should be a boolean
    assert exc.value.args[0] == 'limb should be a boolean value (True/False)'


def test_spherical_limb_invalid5():
    s = SphericalSource()
    with pytest.raises(ValueError) as exc:
        s.limb = 1.  # limb should be a boolean
    assert exc.value.args[0] == 'limb should be a boolean value (True/False)'

# MapSource - specific tests


def test_map_map_none():
    s = MapSource()
    s.map = None


def test_map_map_array():
    s = MapSource()
    s.map = np.ones((2, 3, 4))


def test_map_map_invalid1():
    s = MapSource()
    with pytest.raises(ValueError) as exc:
        s.map = (0., 1., 2., 4.)  # not an array
    assert exc.value.args[0] == 'map should be a Numpy array or an AMRGridView instance'


def test_map_map_invalid2():
    s = MapSource()
    with pytest.raises(ValueError) as exc:
        s.map = [1., 2., 3., 4., 5.]  # not an array
    assert exc.value.args[0] == 'map should be a Numpy array or an AMRGridView instance'


def test_map_map_invalid3():
    s = MapSource()
    with pytest.raises(ValueError) as exc:
        s.map = 'a string'  # not an array
    assert exc.value.args[0] == 'map should be a Numpy array or an AMRGridView instance'


def test_map_write_cartesian():

    # Set up coordinate grid
    g = CartesianGrid([-1., 0., 1.], [-1., 1.], [-1., -0.2, 0.2, 1.])

    v = virtual_file()
    s = MapSource()
    s.luminosity = 1.
    s.map = np.ones((3, 1, 2))
    s.write(v, 'test', g)


def test_map_write_spherical_polar():

    # Set up coordinate grid
    g = SphericalPolarGrid([0., 0.5, 1.], [0., np.pi], [0., 1., 2., 2. * np.pi])

    v = virtual_file()
    s = MapSource()
    s.luminosity = 1.
    s.map = np.ones((3, 1, 2))
    s.write(v, 'test', g)


def test_map_write_cylindrical_polar():

    # Set up coordinate grid
    g = CylindricalPolarGrid([0., 0.5, 1.], [-1., 1.], [0., 1., 2., 2. * np.pi])

    v = virtual_file()
    s = MapSource()
    s.luminosity = 1.
    s.map = np.ones((3, 1, 2))
    s.write(v, 'test', g)


def test_map_write_amr():

    # Set up coordinate grid
    amr = AMRGrid()
    level = amr.add_level()
    grid = level.add_grid()
    grid.xmin, grid.xmax = -1., 1.
    grid.ymin, grid.ymax = -1., 1.
    grid.zmin, grid.zmax = -1., 1.
    grid.nx, grid.ny, grid.nz = 4, 6, 8
    grid.quantities['luminosity'] = np.ones((8, 6, 4))

    v = virtual_file()
    s = MapSource()
    s.luminosity = 1.
    s.map = amr['luminosity']
    s.write(v, 'test', amr)


def test_map_write_octree():

    # Set up coordinate grid
    refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    g = OctreeGrid(0., 0., 0., 10., 10., 10., np.array(refined).astype(bool))

    v = virtual_file()
    s = MapSource()
    s.luminosity = 1.
    s.map = np.ones(len(refined))
    s.write(v, 'test', g)


def test_map_write_invalid1():
    v = virtual_file()
    s = MapSource()
    s.luminosity = 1.
    s.map = None
    with pytest.raises(Exception) as exc:
        s.write(v, 'test', None)
    assert exc.value.args[0] == 'map is not set'


def test_map_write_invalid2():
    v = virtual_file()
    s = MapSource()
    s.luminosity = 1.
    s.map = np.zeros((2, 3, 4))
    with pytest.raises(Exception) as exc:
        s.write(v, 'test', None)
    assert exc.value.args[0] == 'map is zero everywhere'
