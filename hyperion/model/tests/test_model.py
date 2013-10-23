from __future__ import print_function, division

import os
import tempfile
import shutil
from copy import deepcopy

import numpy as np
from astropy.tests.helper import pytest

from .. import Model
from .test_helpers import random_id, get_test_dust, get_realistic_test_dust
from ...grid import CartesianGrid, CylindricalPolarGrid, SphericalPolarGrid, AMRGrid, OctreeGrid
from ...dust import IsotropicDust, SphericalDust

DATA = os.path.join(os.path.dirname(__file__), 'data')


def test_basic(tmpdir):

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.set_n_photons(initial=100, imaging=100)
    m.write(tmpdir.join(random_id()).strpath)


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


class TestAllGridTypes(object):

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
        grid.quantities['density'] = [np.ones((8, 8, 8))]

        refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        self.grid['oct'] = OctreeGrid(0., 0., 0., 10., 10., 10., np.array(refined).astype(bool))

        # Set up initial densities
        self.density = {}
        self.density['car'] = np.array([[[1.]]])
        self.density['cyl'] = np.array([[[1.]]])
        self.density['sph'] = np.array([[[1.]]])
        self.density['amr'] = self.grid['amr']['density'][0]
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
    def test_mismatch_density_energy_2(self, tmpdir, grid_type):
        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust, specific_energy=self.density[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust)
        m.set_n_photons(initial=100, imaging=100)
        with pytest.raises(Exception) as exc:
            m.write(tmpdir.join(random_id()).strpath)
        assert exc.value.args[0] == "Not all dust lists in the grid have the same size"

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_add_density(self, tmpdir, grid_type):
        m = Model()
        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 5000.
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust)
        m.set_n_photons(initial=100, imaging=100)
        m.write(tmpdir.join(random_id()).strpath)
        m.run(tmpdir.join(random_id()).strpath)

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_add_density_from_grid(self, tmpdir, grid_type):
        m = Model()
        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 5000.
        g = deepcopy(self.grid[grid_type])
        if grid_type != 'amr':
            g['density'] = []
            g['density'].append(self.density[grid_type])
        m.set_grid(g)
        print(g['density'])
        m.add_density_grid(g['density'][0], self.dust)
        m.add_density_grid(g['density'][0], self.dust)
        m.set_n_photons(initial=100, imaging=100)
        m.write(tmpdir.join(random_id()).strpath)
        m.run(tmpdir.join(random_id()).strpath)

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_density(self, tmpdir, grid_type):
        m = Model()
        s = m.add_point_source()
        s.luminosity = 1.
        s.temperature = 5000.
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust)
        m.add_density_grid(self.density[grid_type], self.dust, merge_if_possible=True)
        m.set_n_photons(initial=100, imaging=100)
        m.write(tmpdir.join(random_id()).strpath)
        m.run(tmpdir.join(random_id()).strpath)


class TestMerge(object):

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

        self.tmpdir = tempfile.mkdtemp()

        self.dust1_filename = os.path.join(self.tmpdir, random_id())
        self.dust1 = get_test_dust()
        self.dust1.write(self.dust1_filename)

        self.dust2_filename = os.path.join(self.tmpdir, random_id())
        self.dust2 = get_test_dust()
        self.dust2.write(self.dust2_filename)

        self.dust3_filename = os.path.join(self.tmpdir, random_id())
        self.dust3 = IsotropicDust([3.e9, 3.e16], [0.5, 0.5], [1., 0.5])
        self.dust3.emissivities.set_lte(self.dust3.optical_properties, n_temp=10, temp_min=0.1, temp_max=1600.)
        self.dust3.write(self.dust3_filename)

        # The following dust file does not have emissivities and mean
        # opacities since it has never been written to a file
        self.dust4 = get_test_dust()

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_no(self, grid_type):

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust1_filename)
        m.add_density_grid(self.density[grid_type], self.dust2_filename, merge_if_possible=True)
        assert m.grid.n_dust == 2

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_filename_disabled(self, grid_type):

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust1_filename)
        m.add_density_grid(self.density[grid_type], self.dust1_filename, merge_if_possible=False)
        assert m.grid.n_dust == 2

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_filename(self, grid_type):

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust1_filename)
        m.add_density_grid(self.density[grid_type], self.dust1_filename, merge_if_possible=True)
        assert m.grid.n_dust == 1

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_object_identical(self, grid_type):

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust1)
        m.add_density_grid(self.density[grid_type], self.dust1, merge_if_possible=True)
        assert m.grid.n_dust == 1

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_object_samehash(self, grid_type):

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust1)
        m.add_density_grid(self.density[grid_type], self.dust2, merge_if_possible=True)
        assert m.grid.n_dust == 1

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_object_diffhash(self, grid_type):

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust1)
        m.add_density_grid(self.density[grid_type], self.dust3, merge_if_possible=True)
        assert m.grid.n_dust == 2

    @pytest.mark.parametrize(('grid_type'), ['car', 'sph', 'cyl', 'amr', 'oct'])
    def test_merge_object_incomplete(self, grid_type):

        m = Model()
        m.set_grid(self.grid[grid_type])
        m.add_density_grid(self.density[grid_type], self.dust1)
        m.add_density_grid(self.density[grid_type], self.dust4, merge_if_possible=True)
        assert m.grid.n_dust == 2


def test_dust_mix(tmpdir):
    # This is a regression test for a bug which caused the code to crash if
    # isotropic dust and non-isotropic dust were used together.

    iso_dust = get_realistic_test_dust()
    kmh_dust = os.path.join(DATA, 'kmh_lite.hdf5')

    m = Model()

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.

    m.add_density_grid(np.array([[[1.]]]), kmh_dust)
    m.add_density_grid(np.array([[[1.]]]), iso_dust)

    m.set_n_photons(initial=100000, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)
    m.run(tmpdir.join(random_id()).strpath)


def test_dust_changed_nosave(tmpdir):

    kmh_dust = SphericalDust(os.path.join(DATA, 'kmh_lite.hdf5'))
    kmh_dust.set_sublimation_temperature('fast', temperature=1600)

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), kmh_dust)
    m.set_n_photons(initial=1, imaging=1)
    with pytest.raises(ValueError) as exc:
        m.write(tmpdir.join(random_id()).strpath, copy=False)
    assert exc.value.args[0].startswith('Dust properties have been modified since being read in')


def test_dust_changed_save(tmpdir):

    kmh_dust = SphericalDust(os.path.join(DATA, 'kmh_lite.hdf5'))
    kmh_dust.set_sublimation_temperature('fast', temperature=1600)
    kmh_dust.write(tmpdir.join(random_id()).strpath)

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), kmh_dust)
    m.set_n_photons(initial=1, imaging=1)
    m.write(tmpdir.join(random_id()).strpath, copy=False)


def test_model_minimal(tmpdir):

    m = Model()
    m.set_cartesian_grid([-1., 1.],[-1., 1.],[-1., 1.])
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=10)
    m.write(tmpdir.join(random_id()).strpath)
    m.run()
