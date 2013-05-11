from __future__ import print_function, division

from astropy.tests.helper import pytest
import numpy as np

from .. import CartesianGrid, \
               CylindricalPolarGrid, \
               SphericalPolarGrid

GRIDS = [CartesianGrid, CylindricalPolarGrid, SphericalPolarGrid]

WALL = {}
WALL[CartesianGrid] = ['x_wall', 'y_wall', 'z_wall']
WALL[CylindricalPolarGrid] = ['w_wall', 'z_wall', 'p_wall']
WALL[SphericalPolarGrid] = ['r_wall', 't_wall', 'p_wall']


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_list(grid):
    grid([0., 1.], [0., 1.], [0., 1.])


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_tuple(grid):
    grid((0., 1.), (0., 1.), (0., 1.))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_array(grid):
    grid(np.array([0., 1.]),
         np.array([0., 1.]),
         np.array([0., 1.]))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_mixed(grid):
    grid([0., 1.],
         (0., 1.),
         np.array([0., 1.]))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid1(grid):
    with pytest.raises(ValueError) as e:
        grid('hello',  # invalid entry
             (0., 1.),
             (0., 1.))
    assert e.value.args[0] == WALL[grid][0] + ' should be a 1-D sequence'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid2(grid):
    with pytest.raises(ValueError) as e:
        grid((0., 1.),
             1.233,  # invalid entry
             (0., 1.))
    assert e.value.args[0] == WALL[grid][1] + ' should be a 1-D sequence'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid3(grid):
    with pytest.raises(ValueError) as e:
        grid((0., 1.),
             (0., 1.),
             set([1, 2, 3]))  # invalid entry
    assert e.value.args[0] == WALL[grid][2] + ' should be a 1-D sequence'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid4(grid):
    with pytest.raises(ValueError) as e:
        grid([[0., 1.]],  # lists should be 1D
              (0., 1.),
              np.array([0., 1.]))
    assert e.value.args[0] == WALL[grid][0] + ' should be a 1-D sequence'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid5(grid):
    with pytest.raises(ValueError) as e:
        grid([0., 1.],
             ((0., 1.),),  # tuples should be 1D
             np.array([0., 1.]))
    assert e.value.args[0] == WALL[grid][1] + ' should be a 1-D sequence'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid6(grid):
    with pytest.raises(ValueError) as e:
        grid([0., 1.],
             (0., 1.),
             np.array([[0., 1.]]))  # arrays should be 1D
    assert e.value.args[0] == WALL[grid][2] + ' should be a 1-D sequence'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid7(grid):
    with pytest.raises(ValueError) as e:
        grid([-1., -2., 1.],  # should be increasing
             [0., 1.],
             [0., 1.])
    assert e.value.args[0] == WALL[grid][0] + ' should be monotonically increasing'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid8(grid):
    with pytest.raises(ValueError) as e:
        grid([0., 1.],
             [2., -1.],  # should be increasing
             [0., 1.])
    assert e.value.args[0] == WALL[grid][1] + ' should be monotonically increasing'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid9(grid):
    with pytest.raises(ValueError) as e:
        grid([0., 1.],
             [0., 1.],
             [4., -1., 5.])  # should be increasing
    assert e.value.args[0] == WALL[grid][2] + ' should be monotonically increasing'


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_dimension(grid):
    g = grid([0., 1.], [0., 0.5, 1.], [0., 0.25, 0.75, 1.])
    assert g.shape == (3, 2, 1)  # order is reversed


def test_grid_cartesian_ranges():
    g = CartesianGrid([-10., 10.], [-10., 10.], [-10., 10.])


def test_grid_cylindrical_ranges():
    g = CylindricalPolarGrid([-10., 10.], [-10., 10.], [0., 2. * np.pi])


def test_grid_cylindrical_ranges_invalid():
    with pytest.raises(ValueError) as e:
        g = CylindricalPolarGrid([-10., 10.], [-10., 10.], [-10., 10.])
    assert e.value.args[0] == 'p_wall values be in the range [0:2*pi]'


def test_grid_spherical_ranges():
    g = SphericalPolarGrid([-10., 10.], [0., np.pi], [0., 2. * np.pi])


def test_grid_spherical_ranges_invalid1():
    with pytest.raises(ValueError) as e:
        g = SphericalPolarGrid([-10., 10.], [-10., 10.], [0., 2. * np.pi])
    assert e.value.args[0] == 't_wall values be in the range [0:pi]'


def test_grid_spherical_ranges_invalid2():
    with pytest.raises(ValueError) as e:
        g = SphericalPolarGrid([-10., 10.], [0., np.pi], [-10., 10.])
    assert e.value.args[0] == 'p_wall values be in the range [0:2*pi]'
