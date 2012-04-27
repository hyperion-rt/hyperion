from __future__ import print_function, division

import numpy as np
import pytest

from .. import Model
from .test_helpers import random_filename, get_test_dust
from ...grid import CartesianGrid, CylindricalPolarGrid, SphericalPolarGrid, AMRGrid, OctreeGrid


def test_basic():

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.set_n_photons(initial=100, imaging=100)
    m.write(random_filename())


def test_noname_nofilename():
    m = Model()
    with pytest.raises(ValueError) as e:
        m.write()
    assert e.value.args[0] == "filename= has not been specified and model has no name"


def test_nogrid():
    m = Model()
    with pytest.raises(Exception) as e:
        m.write('test')
    assert e.value.args[0] == 'No coordinate grid has been set up'


def test_nophotons():
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    with pytest.raises(Exception) as e:
        m.write('test')
    assert e.value.args[0] == 'Photon numbers not set'


def test_incomplete_photons_1():
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    with pytest.raises(Exception) as e:
        m.set_n_photons(initial=1)
    assert e.value.args[0] == '[n_photons] imaging should bet set'


def test_incomplete_photons_2():
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    with pytest.raises(Exception) as e:
        m.set_n_photons(imaging=1)
    assert e.value.args[0] == '[n_photons] initial should be set since the initial iterations are being computed'


class TestDensitySpecificEnergy(object):

    @classmethod
    def setup_class(self):

        self.grid = {}

        self.grid['car'] = CartesianGrid([-1., 1.], [-2., 2.], [-3., 3.])
        self.grid['cyl'] = CylindricalPolarGrid([0., 1.], [-1., 1.], [0., 2. * np.pi])
        self.grid['sph'] = SphericalPolarGrid([0., 1.], [0., np.pi], [0., 2. * np.pi])

        self.grid['amr'] = AMRGrid()
        level = self.grid['amr'].add_level()
        grid = level.add_grid()
        grid.xmin, grid.xmax = -1., 1.
        grid.ymin, grid.ymax = -1., 1.
        grid.zmin, grid.zmax = -1., 1.
        grid.nx, grid.ny, grid.nz = 8, 8, 8
        grid.quantities['density'] = np.ones((8, 8, 8))

        refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        self.grid['oct'] = OctreeGrid(0., 0., 0., 10., 10., 10., np.array(refined).astype(bool))

        # Set up initial densities
        self.density = {}
        self.density['car'] = np.array([[[1.]]])
        self.density['cyl'] = np.array([[[1.]]])
        self.density['sph'] = np.array([[[1.]]])
        self.density['amr'] = self.grid['amr']['density']
        self.density['oct'] = np.ones(len(refined))

        self.dust = get_test_dust()

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_mismatch_density_energy_1(self, grid_type):
        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust)
        with pytest.raises(Exception) as exc:
            m.add_density_grid(self.density[grid_type], self.dust, specific_energy=self.density[grid_type])
        assert exc.value.args[0] == "Cannot add specific energy as it was not added for previous density arrays"

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_mismatch_density_energy_2(self, grid_type):
        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust, specific_energy=self.density[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust, merge_if_possible=False)
        m.set_n_photons(initial=100, imaging=100)
        with pytest.raises(Exception) as exc:
            m.write(random_filename())
        assert exc.value.args[0] == "Not all dust lists in the grid have the same size"
