from __future__ import print_function, division

import numpy as np
import pytest

from .. import Model
from .test_helpers import random_filename, get_test_dust
from ...grid import CartesianGrid, CylindricalPolarGrid, SphericalPolarGrid, AMRGrid, OctreeGrid
from ...dust import IsotropicDust


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
        m.add_density_grid(self.density[grid_type], self.dust)
        m.set_n_photons(initial=100, imaging=100)
        with pytest.raises(Exception) as exc:
            m.write(random_filename())
        assert exc.value.args[0] == "Not all dust lists in the grid have the same size"


class TestMerge(object):

    @classmethod
    def setup_class(self):

        self.dust1_filename = random_filename()
        self.dust1 = get_test_dust()
        self.dust1.write(self.dust1_filename)

        self.dust2_filename = random_filename()
        self.dust2 = get_test_dust()
        self.dust2.write(self.dust2_filename)

        self.dust3_filename = random_filename()
        self.dust3 = IsotropicDust([3.e9, 3.e16], [0.5, 0.5], [1., 0.5])
        self.dust3.emissivities.set_lte(self.dust3.optical_properties, n_temp=10, temp_min=0.1, temp_max=1600.)
        self.dust3.write(self.dust3_filename)

        # The following dust file does not have emissivities and mean
        # opacities since it has never been written to a file
        self.dust4 = get_test_dust()


    def test_merge_no(self):

        m = Model()
        m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
        m.set_n_photons(initial=100, imaging=100)
        m.add_density_grid(np.array([[[1.]]]), self.dust1_filename)
        m.add_density_grid(np.array([[[1.]]]), self.dust2_filename, merge_if_possible=True)
        assert m.grid.n_dust == 2

    def test_merge_filename(self):

        m = Model()
        m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
        m.set_n_photons(initial=100, imaging=100)
        m.add_density_grid(np.array([[[1.]]]), self.dust1_filename)
        m.add_density_grid(np.array([[[1.]]]), self.dust1_filename, merge_if_possible=True)
        assert m.grid.n_dust == 1

    def test_merge_object_identical(self):

        m = Model()
        m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
        m.set_n_photons(initial=100, imaging=100)
        m.add_density_grid(np.array([[[1.]]]), self.dust1)
        m.add_density_grid(np.array([[[1.]]]), self.dust1, merge_if_possible=True)
        assert m.grid.n_dust == 1

    def test_merge_object_samehash(self):

        m = Model()
        m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
        m.set_n_photons(initial=100, imaging=100)
        m.add_density_grid(np.array([[[1.]]]), self.dust1)
        m.add_density_grid(np.array([[[1.]]]), self.dust2, merge_if_possible=True)
        assert m.grid.n_dust == 1

    def test_merge_object_diffhash(self):

        m = Model()
        m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
        m.set_n_photons(initial=100, imaging=100)
        m.add_density_grid(np.array([[[1.]]]), self.dust1)
        m.add_density_grid(np.array([[[1.]]]), self.dust3, merge_if_possible=True)
        assert m.grid.n_dust == 2

    def test_merge_object_incomplete(self):

        m = Model()
        m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
        m.set_n_photons(initial=100, imaging=100)
        m.add_density_grid(np.array([[[1.]]]), self.dust1)
        m.add_density_grid(np.array([[[1.]]]), self.dust4, merge_if_possible=True)
        assert m.grid.n_dust == 2
