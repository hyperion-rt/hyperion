import os

import h5py
import numpy as np

from ..sph import construct_octree

DATA = os.path.join(os.path.dirname(__file__), 'data')


def test_construct_octree():

    np.random.seed(0)

    N = 5000
    px = np.random.uniform(-10., 10., N)
    py = np.random.uniform(-10., 10., N)
    pz = np.random.uniform(-10., 10., N)
    mass = np.random.uniform(0., 1., N)
    sigma = np.random.uniform(0., 0.1, N)

    o = construct_octree(0.1, 0.2, 0.3, 6., 5., 4.,
                         px, py, pz,
                         sigma, mass, n_levels=10)

    # The following lines can be used to write out the reference file if the
    # SPH gridding code is updated.
    # f = h5py.File('reference_octree.hdf5', 'w')
    # o.write(f)
    # f.close()

    from hyperion.grid import OctreeGrid
    f = h5py.File(os.path.join(DATA, 'reference_octree.hdf5'), 'r')
    o_ref = OctreeGrid()
    o_ref.read(f)
    f.close()

    assert np.all(o_ref.refined == o.refined)
    assert np.all(o_ref['density'][0].array == o['density'][0].array)
