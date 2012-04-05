import pytest
import numpy as np

from hyperion.grid import CartesianGrid, \
                          CylindricalPolarGrid, \
                          SphericalPolarGrid

GRIDS = [CartesianGrid, CylindricalPolarGrid, SphericalPolarGrid]


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_list(grid):
    grid([-1., 1.], [-1., 1.], [-1., 1.])


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_tuple(grid):
    grid((-1., 1.), (-1., 1.), (-1., 1.))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_array(grid):
    grid(np.array([-1., 1.]),
         np.array([-1., 1.]),
         np.array([-1., 1.]))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_mixed(grid):
    grid([-1., 1.],
         (-1., 1.),
         np.array([-1., 1.]))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid1(grid):
    with pytest.raises(ValueError):
        grid('hello',  # invalid entry
             (-1., 1.),
             (-1., 1.))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid2(grid):
    with pytest.raises(ValueError):
        grid((-1., 1.),
             1.233,  # invalid entry
             (-1., 1.))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid3(grid):
    with pytest.raises(ValueError):
        grid((-1., 1.),
             (-1., 1.),
             set([1, 2, 3]))  # invalid entry


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid4(grid):
    with pytest.raises(ValueError):
        grid([[-1., 1.]],  # lists should be 1D
              (-1., 1.),
              np.array([-1., 1.]))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid5(grid):
    with pytest.raises(ValueError):
        grid([-1., 1.],
             ((-1., 1.),),  # tuples should be 1D
             np.array([-1., 1.]))


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid6(grid):
    with pytest.raises(ValueError):
        grid([-1., 1.],
             (-1., 1.),
             np.array([[-1., 1.]]))  # arrays should be 1D


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid7(grid):
    with pytest.raises(ValueError):
        grid([-1., -2., 1.],  # should be increasing
             [-1., 1.],
             [-1., 1.])


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid8(grid):
    with pytest.raises(ValueError):
        grid([-1., 1.],
             [2., -1.],  # should be increasing
             [-1., 1.])


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_invalid9(grid):
    with pytest.raises(ValueError):
        grid([-1., 1.],
             [-1., 1.],
             [4., -1., 5.])  # should be increasing


@pytest.mark.parametrize(('grid'), GRIDS)
def test_grid_dimension(grid):
    g = grid([-1., 1.], [-1., 0., 1.], [-1., -0.2, 0.2, 1.])
    assert g.shape == (3, 2, 1)  # order is reversed
