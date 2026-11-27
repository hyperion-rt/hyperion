
import os
import random

import h5py
import numpy as np
from astropy.table import Table

from .. import voronoi_helpers as vh
from ..voronoi_grid import VoronoiGrid


DATA = os.path.join(os.path.dirname(__file__), 'data')

# Helper functions.


# Check that all non-wall neighbours in voro are also in qhull. qhull
# could have extra neighbours due to the clipping of the domain
# by voro++.
def _neighbours_check(voro, qhull):
    idx = voro.st[1]
    vneighs = voro.st[0]
    for i in range(len(idx) - 1):
        cur_neighs = filter(lambda n: n >= 0, vneighs[idx[i]:idx[i + 1]])
        if not all([_ in qhull['neighbours'][i] for _ in cur_neighs]):
            return False
    return True

# Check that the volumes are positive and amount to the total domain volume.


def _volume_check(voro):
    c1 = all(v > 1E-13 for v in voro['volume'])
    c2 = abs(sum(voro['volume']) - 1.) < 1E-13
    return (c1 and c2)

# Check that the bounding boxes are consistent with the sites' positions.


def _bb_check(voro):
    retval = []
    for r in voro:
        c, m, M = r['coordinates'], r['bb_min'], r['bb_max']
        retval.append(all(c[i] >= m[i] and c[i] <= M[i] for i in range(0, 3)))
    return all(retval)

# Test for internal consistency and consistency wrt
# the qhull version of the code.


def test_consistency():
    # Three pickled neighbours tables generated with qhull.
    l = ['qhull_00.hdf5', 'qhull_01.hdf5', 'qhull_02.hdf5']
    for s in l:
        # The three qhull tables have been generated with seeds (0,1,2)
        # and with (10,100,1000) sites respectively.
        idx = l.index(s)
        random.seed(idx)
        sites_arr = np.array([(random.uniform(0, 1), random.uniform(
            0, 1), random.uniform(0, 1)) for i in range(0, 10 ** (idx + 1))])
        vg = vh.voronoi_grid(sites_arr, np.array([[0, 1.], [0, 1], [0, 1]]))
        voro = vg.neighbours_table
        qhull = Table.read(os.path.join(DATA, s), path='data')
        assert _neighbours_check(vg, qhull)
        assert _volume_check(voro)
        assert _bb_check(voro)


def test_evaluate_function_average():

    np.random.seed(12345)

    N = 10
    x = np.random.random(N)
    y = np.random.random(N)
    z = np.random.random(N)

    def gaussian(x, y, z):
        return np.exp(-(x ** 2 + y ** 2 + z ** 2))

    g = VoronoiGrid(x, y, z)
    dens = g.evaluate_function_average(gaussian, n_samples=10000000, min_cell_samples=5)

    r = np.sqrt(x * x + y * y + z * z)

    np.testing.assert_allclose(r, [1.355005, 1.240542, 0.502539, 0.398975, 0.777704,
                                   1.14804, 1.604416, 1.282766, 1.454925, 0.989993], rtol=1e-5)

    np.testing.assert_allclose(dens, [0.202441, 0.29333, 0.522314, 0.767736, 0.451614,
                                      0.294537, 0.121397, 0.224867, 0.124017, 0.373834], rtol=1e-2)


def test_io(tmpdir):

    file_1 = tmpdir.join('test1.hdf5').strpath
    file_2 = tmpdir.join('test2.hdf5').strpath

    g1 = VoronoiGrid([1, 2, 3], [3, 4, 5], [2, 3, 4])
    f1 = h5py.File(file_1, 'w')
    g1.write(f1)
    f1.close()

    f2 = h5py.File(file_1, 'r')
    f3 = h5py.File(file_2, 'w')
    g2 = VoronoiGrid()
    g2.read(f2)
    g2.write(f3)

    np.testing.assert_allclose(g1._sparse_neighbors[0], g2._sparse_neighbors[0])
    np.testing.assert_allclose(g1._sparse_neighbors[1], g2._sparse_neighbors[1])



def test_init_sparse(tmpdir):

    # Regression test for a bug that caused VoronoiGrid._sparse_neighbors to not be defined
    # when writing out a Voronoi grid.

    file_1 = tmpdir.join('test1.hdf5').strpath
    file_2 = tmpdir.join('test2.hdf5').strpath

    g1 = VoronoiGrid([1,2,3],[3,4,5],[2,3,4])
    f1 = h5py.File(file_1, 'w')
    g1.write(f1)
    f1.close()

    f2 = h5py.File(file_1, 'r')
    g2 = VoronoiGrid()
    g2.read(f2)

    f3 = h5py.File(file_2, 'w')
    g3 = VoronoiGrid(g2)
    g3.write(f3)



def test_setitem_from_view():
    # Regression test for VoronoiGrid.__setitem__, which used to check a
    # nonexistent 'refined' attribute when assigning a view to an empty grid
    g1 = VoronoiGrid([1., 2., 3.], [3., 4., 5.], [2., 3., 4.])
    g1.quantities['density'] = np.array([1., 2., 3.])
    g2 = VoronoiGrid()
    g2['density'] = g1['density']
    np.testing.assert_allclose(g2['density'].array, g1['density'].array)
    assert g2.xmin == g1.xmin and g2.xmax == g1.xmax
    assert g2.ymin == g1.ymin and g2.ymax == g1.ymax
    assert g2.zmin == g1.zmin and g2.zmax == g1.zmax


def test_evaluate_function_average_constant():
    # For a constant function every (non-empty) cell average must equal that
    # constant. The old code divided the last cell's sum by a decremented
    # sample count, inflating that one cell's average away from the constant.
    np.random.seed(9876)
    N = 12
    x = np.random.random(N)
    y = np.random.random(N)
    z = np.random.random(N)

    g = VoronoiGrid(x, y, z)
    avg = g.evaluate_function_average(lambda x, y, z: 3.5 * np.ones_like(x),
                                      n_samples=2000000, min_cell_samples=5)
    np.testing.assert_allclose(avg, 3.5, rtol=1e-6)


def _sites_in_sphere(n, x0, y0, z0, r0, seed=42):
    # Random sites uniformly distributed inside a sphere
    rng = np.random.RandomState(seed)
    r = r0 * rng.uniform(0., 1., n) ** (1. / 3.)
    costh = rng.uniform(-1., 1., n)
    phi = rng.uniform(0., 2. * np.pi, n)
    sinth = np.sqrt(1. - costh ** 2)
    x = x0 + r * sinth * np.cos(phi)
    y = y0 + r * sinth * np.sin(phi)
    z = z0 + r * costh
    return x, y, z


def test_voronoi_sphere_wall_volumes():
    # With a spherical wall, the cells at the edge of the site distribution
    # are bounded by the wall instead of extending to the domain box, so the
    # total volume should be close to the volume of the sphere. Since voro++
    # cuts each cell with the plane tangent to the sphere, the sphere is
    # circumscribed by planar facets and the total volume converges to the
    # sphere volume from above as the number of sites increases.
    x, y, z = _sites_in_sphere(10000, 1., 2., 3., 5.)
    sites = np.array([x, y, z]).transpose()
    domain = np.array([[-6., 8.], [-5., 9.], [-4., 10.]])

    g_nowall = vh.voronoi_grid(sites, domain)
    volume_nowall = np.sum(g_nowall.neighbours_table['volume'])
    np.testing.assert_allclose(volume_nowall, 14. ** 3)

    g_wall = vh.voronoi_grid(sites, domain, wall='sphere', wall_args=(1., 2., 3., 5.))
    volume_wall = np.sum(g_wall.neighbours_table['volume'])
    volume_sphere = 4. / 3. * np.pi * 5. ** 3
    assert volume_sphere < volume_wall < 1.05 * volume_sphere

    # The wall faces should appear in the neighbour lists with the id -7
    assert np.min(g_wall.st[0]) == -7


def test_voronoi_sphere_wall_invalid():
    import pytest
    x, y, z = _sites_in_sphere(100, 0., 0., 0., 5.)
    sites = np.array([x, y, z]).transpose()
    domain = np.array([[-6., 6.], [-6., 6.], [-6., 6.]])
    with pytest.raises(ValueError):
        vh.voronoi_grid(sites, domain, wall='cube', wall_args=(0., 0., 0., 5.))
    with pytest.raises(ValueError):
        vh.voronoi_grid(sites, domain, wall='sphere', wall_args=(0., 0., 5.))
    with pytest.raises(ValueError):
        vh.voronoi_grid(sites, domain, wall='sphere', wall_args=(0., 0., 0., -5.))
    with pytest.raises(ValueError):
        # sites outside the wall
        vh.voronoi_grid(sites, domain, wall='sphere', wall_args=(0., 0., 0., 1.))


def test_voronoi_grid_sphere_wall_roundtrip(tmpdir):
    # The wall parameters should survive a write/read round-trip, including
    # the geometry hash check on read.
    x, y, z = _sites_in_sphere(1000, 0., 0., 0., 5.)
    g = VoronoiGrid(x, y, z,
                    xmin=-6., xmax=6., ymin=-6., ymax=6., zmin=-6., zmax=6.,
                    wall='sphere', wall_args=(0., 0., 0., 5.))
    volume_sphere = 4. / 3. * np.pi * 5. ** 3
    assert np.sum(g.volumes) < 1.2 * volume_sphere

    filename = os.path.join(str(tmpdir), 'wall_grid.h5')
    with h5py.File(filename, 'w') as f:
        g.write(f)
    with h5py.File(filename, 'r') as f:
        g2 = VoronoiGrid()
        g2.read(f)
    assert g2._wall == 'sphere'
    assert g2._wall_args == (0., 0., 0., 5.)
    np.testing.assert_allclose(np.sum(g2.volumes), np.sum(g.volumes))
