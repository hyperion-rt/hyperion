import os
import sys
from copy import deepcopy
import h5py
import numpy as np
import pytest

from ...util.functions import random_id
from ...model import Model
from ...model.tests.test_helpers import get_realistic_test_dust

from .. import (CartesianGrid,
                SphericalPolarGrid,
                CylindricalPolarGrid,
                AMRGrid,
                OctreeGrid)

try:
    import yt
except:
    YT_VERSION = None
else:
    YT_VERSION = 3

DATA = os.path.join(os.path.dirname(__file__), 'data')

ALL_GRID_TYPES = ['car', 'amr', 'oct']


@pytest.mark.skipif("YT_VERSION is None")
class TestToYt(object):

    def setup_method(self, method):

        # Set up grids
        self.grid = {}

        self.grid['car'] = CartesianGrid([-1., 1.],
                                         [-2., 0., 2.],
                                         [-3., -1., 1., 3.])
        self.grid['cyl'] = CylindricalPolarGrid([0., 1.],
                                                [-1., 0., 1.],
                                                [0., 0.75 * np.pi, 1.25 * np.pi, 2. * np.pi])
        self.grid['sph'] = SphericalPolarGrid([0., 0.5, 1.],
                                              [0., np.pi],
                                              [0., 0.75 * np.pi, 1.25 * np.pi, 2. * np.pi])

        self.grid['amr'] = AMRGrid()
        level = self.grid['amr'].add_level()
        grid = level.add_grid()
        grid.xmin, grid.xmax = -1., 1.
        grid.ymin, grid.ymax = -1., 1.
        grid.zmin, grid.zmax = -1., 1.
        grid.nx, grid.ny, grid.nz = 8, 6, 4

        refined = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        self.grid['oct'] = OctreeGrid(0., 0., 0., 10., 10., 10.,
                                      np.array(refined).astype(bool))

        # Set up initial densities
        self.density = {}
        self.density['car'] = np.ones((3, 2, 1))
        self.density['cyl'] = np.ones((3, 2, 1))
        self.density['sph'] = np.ones((3, 1, 2))
        amr_q = deepcopy(self.grid['amr'])
        amr_q.levels[0].grids[0].quantities['density'] = np.ones((4, 6, 8))
        self.density['amr'] = amr_q['density']
        self.density['oct'] = np.ones(len(refined))

    @pytest.mark.parametrize(('grid_type'), ALL_GRID_TYPES)
    def test_to_yt(self, tmpdir, grid_type):

        from yt import ProjectionPlot

        g = self.grid[grid_type]
        g['density'] = []
        g['density'].append(self.density[grid_type])
        pf = g.to_yt()

        p = ProjectionPlot(pf, 'x', ["density"], center='c', origin='native')
        p.save(tmpdir.join('test.png').strpath)


@pytest.mark.requires_hyperion_binaries
@pytest.mark.skipif("YT_VERSION is None")
def test_from_yt(tmpdir):

    from yt import load

    ds = load(os.path.join(DATA, 'DD0010', 'moving7_0010'))

    def _dust_density(field, data):
        return data["density"].in_units('g/cm**3') * 0.01

    ds.add_field(("gas", "dust_density"), function=_dust_density, units='g/cm**3', sampling_type='cell')

    amr = AMRGrid.from_yt(ds, quantity_mapping={'density': ('gas', 'dust_density')})

    m = Model()

    m.set_amr_grid(amr)

    m.add_density_grid(amr['density'], get_realistic_test_dust())

    s = m.add_point_source()
    s.luminosity = 1000
    s.temperature = 1000

    m.set_n_initial_iterations(3)

    m.set_n_photons(initial=1e5, imaging=0)

    m.set_propagation_check_frequency(1)

    m.set_copy_input(False)

    input_file = tmpdir.join('test.rtin').strpath
    output_file = tmpdir.join('test.rtout').strpath

    m.write(input_file)
    m.run(output_file)


@pytest.mark.skipif("YT_VERSION is None")
def test_axis_ordering_cartesian():

    # Regression test for axis ordering

    from .yt_compat import get_frb

    x = np.linspace(-1, 1, 9)
    y = np.linspace(-2, 2, 17)
    z = np.linspace(-3, 3, 33)

    density = np.arange(32)[:, None, None] * np.ones((32, 16, 8))

    g = CartesianGrid(x, y, z)
    g['density'] = []
    g['density'].append(density)

    from yt import ProjectionPlot, SlicePlot

    pf = g.to_yt()

    for iz, z in enumerate(g.z):
        prj = SlicePlot(pf, 'z', ['density'], center=[0.0, 0.0, z])
        np.testing.assert_allclose(get_frb(prj, 'density').min(), iz)
        np.testing.assert_allclose(get_frb(prj, 'density').max(), iz)


@pytest.mark.skipif("YT_VERSION is None")
def test_axis_ordering_amr():

    # Regression test for axis ordering

    from .yt_compat import get_frb

    g = AMRGrid()

    level = g.add_level()

    grid = level.add_grid()
    grid.xmin, grid.xmax = -1, 1
    grid.ymin, grid.ymax = -2, 2
    grid.zmin, grid.zmax = -3, 3
    grid.nx, grid.ny, grid.nz = 8, 16, 32

    grid.quantities['density'] = []
    grid.quantities['density'].append(np.arange(grid.nz)[:, None, None] * np.ones((grid.nz, grid.ny, grid.nx)))

    from yt import ProjectionPlot, SlicePlot

    pf = g.to_yt()

    zw = np.linspace(-3, 3, 33)
    zcen = 0.5 * (zw[1:] + zw[:-1])

    for iz, z in enumerate(zcen):
        prj = SlicePlot(pf, 'z', ['density'], center=[0.0, 0.0, z])
        np.testing.assert_allclose(get_frb(prj, 'density').min(), iz)
        np.testing.assert_allclose(get_frb(prj, 'density').max(), iz)


@pytest.mark.skipif("YT_VERSION is None")
def test_octree_to_yt_positions():

    # Regression test for a bug that caused octree grids to fail to load in
    # yt 4.1 and later, which also checks that each cell in yt has the same
    # position and value as the corresponding cell in the octree

    refined = np.array([1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]).astype(bool)

    g = OctreeGrid(1., 2., 3., 4., 6., 8., refined)
    density = np.arange(len(refined)) + 1.
    g['density'] = []
    g['density'].append(density)

    ds = g.to_yt()

    np.testing.assert_allclose(ds.domain_left_edge.ndarray_view(), [-3., -4., -5.])
    np.testing.assert_allclose(ds.domain_right_edge.ndarray_view(), [5., 8., 11.])

    ad = ds.all_data()
    x, y, z = [ad['index', axis].ndarray_view() for axis in 'xyz']
    values = ad['gas', 'density'].ndarray_view()

    xmin, xmax, ymin, ymax, zmin, zmax = g.limits
    leaves = ~refined
    expected = {(0.5 * (xmin[i] + xmax[i]), 0.5 * (ymin[i] + ymax[i]), 0.5 * (zmin[i] + zmax[i])): density[i]
                for i in np.nonzero(leaves)[0]}

    assert len(values) == np.sum(leaves)
    for xc, yc, zc, value in zip(x, y, z, values):
        assert expected[(xc, yc, zc)] == value
