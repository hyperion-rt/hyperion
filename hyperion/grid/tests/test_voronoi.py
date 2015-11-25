from __future__ import print_function, division

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

