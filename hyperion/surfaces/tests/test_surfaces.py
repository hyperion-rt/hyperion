from __future__ import print_function, division

import string
import random

import h5py
import pytest
import numpy as np

from ..surfaces import Surface, SphericalSurface

from ..surface_scattering_properties import SurfaceScatteringProperties

ALL_SURFACES = [Surface, SphericalSurface]


def random_id(length=32):
    return ''.join(random.sample(string.ascii_letters + string.digits, length))


def virtual_file():
    return h5py.File(random_id(), driver='core', backing_store=False)

@pytest.mark.parametrize(('surface_type'), ALL_SURFACES)
def test_init(surface_type):
    surface_type()

# SURFACE PROPERTIES

@pytest.mark.parametrize(('surface_type'), ALL_SURFACES)
def test_surface_properties_class(surface_type):
    s = surface_type()
    s.surface_properties = SurfaceScatteringProperties()


@pytest.mark.parametrize(('surface_type'), ALL_SURFACES)
def test_surface_properties_string(surface_type):
    s = surface_type()
    s.surface_properties = 'filename'


@pytest.mark.parametrize(('surface_type'), ALL_SURFACES)
def test_surface_properties_invalid(surface_type):
    s = surface_type()
    with pytest.raises(TypeError) as exc:
        s.surface_properties = 3.14
    assert exc.value.args[0] == 'surface_properties should be a string or a SurfaceScatteringProperties instance'


# # POSITION
#
SURFACES_POSITION = [SphericalSurface]


def test_position_tests_complete():
    expected = [x for x in ALL_SURFACES if hasattr(x, 'position')]
    extra = list(set(SURFACES_POSITION) - set(expected))
    assert extra == []


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_none(surface_type):
    s = surface_type()
    s.position = None


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_tuple(surface_type):
    s = surface_type()
    s.position = (0., 1., 2.)


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_tuple_invalid(surface_type):
    s = surface_type()
    with pytest.raises(ValueError) as exc:
        s.position = (0., 1., 2., 4.)  # too many elements
    assert exc.value.args[0] == 'position should be a sequence of 3 values'


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_list(surface_type):
    s = surface_type()
    s.position = [1., 2., 3.]


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_list_invalid(surface_type):
    s = surface_type()
    with pytest.raises(ValueError) as exc:
        s.position = [1., 2.]  # too few elements
    assert exc.value.args[0] == 'position should be a sequence of 3 values'


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_numpy(surface_type):
    s = surface_type()
    s.position = np.array([2., 3., 4.])


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_numpy_invalid1(surface_type):
    s = surface_type()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([2.])  # too few elements
    assert exc.value.args[0] == 'position should be a sequence of 3 values'


@pytest.mark.parametrize(('surface_type'), SURFACES_POSITION)
def test_position_numpy_invalid2(surface_type):
    s = surface_type()
    with pytest.raises(ValueError) as exc:
        s.position = np.array([[1., 2., 3.]])  # wrong dimensionality
    assert exc.value.args[0] == 'position should be a 1-D sequence'

# RADIUS

SURFACES_RADIUS = [SphericalSurface]


def test_radius_tests_complete():
    expected = [x for x in ALL_SURFACES if hasattr(x, 'radius')]
    extra = list(set(SURFACES_RADIUS) - set(expected))
    assert extra == []


@pytest.mark.parametrize(('surface_type'), SURFACES_RADIUS)
def test_radius_none(surface_type):
    s = surface_type()
    s.radius = None


@pytest.mark.parametrize(('surface_type'), SURFACES_RADIUS)
def test_radius_float(surface_type):
    s = surface_type()
    s.radius = 1.e10


@pytest.mark.parametrize(('surface_type'), SURFACES_RADIUS)
def test_radius_invalid2(surface_type):
    s = surface_type()
    with pytest.raises(ValueError) as exc:
        s.radius = np.array([1, 2, 3])  # radius should be a scalar
    assert exc.value.args[0] == 'radius should be a scalar value'


@pytest.mark.parametrize(('surface_type'), SURFACES_RADIUS)
def test_radius_invalid3(surface_type):
    s = surface_type()
    with pytest.raises(ValueError) as exc:
        s.radius = 'invalid'  # radius should be a number
    assert exc.value.args[0] == 'radius should be a numerical value'


@pytest.mark.parametrize(('surface_type'), SURFACES_RADIUS)
def test_radius_invalid4(surface_type):
    s = surface_type()
    with pytest.raises(ValueError) as exc:
        s.radius = -1.  # radius should be positive
    assert exc.value.args[0] == 'radius should be positive'

example_properties = SurfaceScatteringProperties()
example_properties.nu = [0.1, 0.5, 0.8]
example_properties.mu0 = [0., 1.]
example_properties.mu = [0., 0.33, 0.66, 1.0]
example_properties.psi = [0., 1., 2., 3., 4.]
example_properties.albedo = [0.2, 0.3, 0.4]
example_properties.brdf = np.arange(120).reshape(3, 2, 4, 5)

REQUIRED = {}

REQUIRED[Surface] = {}
REQUIRED[Surface]['surface_properties'] = example_properties

REQUIRED[SphericalSurface] = {}
REQUIRED[SphericalSurface]['radius'] = 1.
REQUIRED[SphericalSurface]['surface_properties'] = example_properties

# Note that position is not required since it defaults to the origin.
# limb is also not required, since it defaults to False.

# Test that no errors are raised if all attributes are present


@pytest.mark.parametrize(('surface_type'), REQUIRED)
def test_all(surface_type):
    v = virtual_file()
    s = surface_type()
    for attribute in REQUIRED[surface_type]:
        setattr(s, attribute, REQUIRED[surface_type][attribute])
    if surface_type == Surface:
        s.write(v)
    else:
        s.write(v, 'test')

# Test that an error is raised if one attribute is missing

combinations = [(s, a) for s in REQUIRED for a in REQUIRED[s]]


@pytest.mark.parametrize(('surface_type', 'missing'), combinations)
def test_missing(surface_type, missing):
    v = virtual_file()
    s = surface_type()
    for attribute in REQUIRED[surface_type]:
        if attribute == missing:
            continue
        setattr(s, attribute, REQUIRED[surface_type][attribute])
    with pytest.raises(ValueError) as exc:
        if surface_type == Surface:
            s.write(v)
        else:
            s.write(v, 'test')
    assert exc.value.args[0] == '{0:s} is not set'.format(missing)
