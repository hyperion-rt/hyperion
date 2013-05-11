from copy import deepcopy

import h5py
import numpy as np
from astropy.tests.helper import pytest

from ...util.functions import random_id
from .. import CartesianGrid, \
               CylindricalPolarGrid, \
               SphericalPolarGrid, \
               AMRGrid, \
               OctreeGrid

ALL_GRID_TYPES = ['car', 'sph', 'cyl', 'amr', 'oct']


def exc_msg(exc):
    if isinstance(exc.value, basestring):
        return exc.value
    elif type(exc.value) is tuple:
        return exc.value[0]
    else:
        return exc.value.args[0]


class TestView(object):

    def setup_method(self, method):

        # Set up grids
        self.grid = {}

        self.grid['car'] = CartesianGrid([-1., 1.],
                                         [-2., 2.],
                                         [-3., 3.])
        self.grid['cyl'] = CylindricalPolarGrid([0., 1.],
                                                [-1., 1.],
                                                [0., 2. * np.pi])
        self.grid['sph'] = SphericalPolarGrid([0., 1.],
                                              [0., np.pi],
                                              [0., 2. * np.pi])

        self.grid['amr'] = AMRGrid()
        level = self.grid['amr'].add_level()
        grid = level.add_grid()
        grid.xmin, grid.xmax = -1., 1.
        grid.ymin, grid.ymax = -1., 1.
        grid.zmin, grid.zmax = -1., 1.
        grid.nx, grid.ny, grid.nz = 8, 8, 8

        refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        self.grid['oct'] = OctreeGrid(0., 0., 0., 10., 10., 10.,
                                      np.array(refined).astype(bool))

        # Set up empty grid class
        self.grid_empty = {}
        self.grid_empty['car'] = CartesianGrid
        self.grid_empty['cyl'] = CylindricalPolarGrid
        self.grid_empty['sph'] = SphericalPolarGrid
        self.grid_empty['amr'] = AMRGrid
        self.grid_empty['oct'] = OctreeGrid

        # Set up initial densities
        self.density = {}
        self.density['car'] = np.array([[[1.]]])
        self.density['cyl'] = np.array([[[1.]]])
        self.density['sph'] = np.array([[[1.]]])
        amr_q = deepcopy(self.grid['amr'])
        amr_q.levels[0].grids[0].quantities['density'] = np.ones((8, 8, 8))
        self.density['amr'] = amr_q['density']
        self.density['oct'] = np.ones(len(refined))

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_write_read_empty(self, grid_type):
        g = self.grid[grid_type]
        f = h5py.File(random_id(), driver='core', backing_store=False)
        g.write(f)
        h = self.grid_empty[grid_type]()
        h.read(f)
        f.close()
        assert h.n_dust is None

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_write_read_single(self, grid_type):
        g = self.grid[grid_type]
        f = h5py.File(random_id(), driver='core', backing_store=False)
        g['density'] = []
        g['density'].append(self.density[grid_type])
        g.write(f)
        h = self.grid_empty[grid_type]()
        h.read(f)
        f.close()
        assert h.n_dust == 1
        if grid_type == 'amr':
            assert type(h.levels[0].grids[0].quantities['density']) is list
        else:
            assert type(h.quantities['density']) is list

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_write_read_double(self, grid_type):
        g = self.grid[grid_type]
        f = h5py.File(random_id(), driver='core', backing_store=False)
        g['density'] = []
        g['density'].append(self.density[grid_type])
        g['density'].append(self.density[grid_type])
        g.write(f)
        h = self.grid_empty[grid_type]()
        h.read(f)
        f.close()
        assert h.n_dust == 2
        if grid_type == 'amr':
            assert type(h.levels[0].grids[0].quantities['density']) is list
        else:
            assert type(h.quantities['density']) is list

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_write_read_double_multiple(self, grid_type):
        g = self.grid[grid_type]
        f = h5py.File(random_id(), driver='core', backing_store=False)
        g['density'] = []
        g['density'].append(self.density[grid_type])
        g['density'].append(self.density[grid_type])
        g['energy'] = []
        g['energy'].append(self.density[grid_type])
        g['energy'].append(self.density[grid_type])
        g.write(f)
        h = self.grid_empty[grid_type]()
        h.read(f)
        f.close()
        assert h.n_dust == 2
        if grid_type == 'amr':
            assert type(h.levels[0].grids[0].quantities['density']) is list
            assert type(h.levels[0].grids[0].quantities['energy']) is list
        else:
            assert type(h.quantities['density']) is list
            assert type(h.quantities['energy']) is list

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_write_read_type_mismatch(self, grid_type):
        g = self.grid[grid_type]
        f = h5py.File(random_id(), driver='core', backing_store=False)
        g['density'] = []
        g['density'].append(self.density[grid_type])
        g.write(f)
        f['Geometry'].attrs['grid_type'] = 'invalid'.encode('utf-8')
        h = self.grid_empty[grid_type]()
        with pytest.raises(Exception) as exc:
            h.read(f)
        if grid_type == 'car':
            assert exc.value.args[0] == "Grid is not cartesian"
        elif grid_type == 'cyl':
            assert exc.value.args[0] == "Grid is not cylindrical polar"
        elif grid_type == 'sph':
            assert exc.value.args[0] == "Grid is not spherical polar"
        elif grid_type == 'amr':
            assert exc.value.args[0] == "Grid is not an AMR grid"
        elif grid_type == 'oct':
            assert exc.value.args[0] == "Grid is not an octree"

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_write_read_hash_mismatch(self, grid_type):
        g = self.grid[grid_type]
        f = h5py.File(random_id(), driver='core', backing_store=False)
        g['density'] = []
        g['density'].append(self.density[grid_type])
        g.write(f)
        f['Geometry'].attrs['geometry'] = 'a4e2805a72dfcf01b2fd94da0be32511'.encode('utf-8')
        h = self.grid_empty[grid_type]()
        with pytest.raises(Exception) as exc:
            h.read(f)
        assert exc.value.args[0] == "Calculated geometry hash does not " \
                                    "match hash in file"

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_write_read_groups_exist(self, grid_type):
        g = self.grid[grid_type]
        f = h5py.File(random_id(), driver='core', backing_store=False)
        f.create_group('Geometry')
        f.create_group('Quantities')
        g['density'] = []
        g['density'].append(self.density[grid_type])
        g.write(f)
        h = self.grid_empty[grid_type]()
        h.read(f)
        assert h.n_dust == 1
