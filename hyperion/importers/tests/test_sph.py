import os

import h5py
import numpy as np
from numpy.testing import assert_allclose

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
    assert_allclose(o_ref['density'][0].array, o['density'][0].array)

    # Find limits
    xmin, xmax, ymin, ymax, zmin, zmax = o.limits

    # Find volumes
    volumes = o.volumes

    # Check that volume computed from limits matches
    assert_allclose(volumes, (xmax - xmin) * (ymax - ymin) * (zmax - zmin))

    # Now we check that the volumes computed by Hyperion are equivalent to ones
    # computed using a recursive (and slow) technique. This used to be what
    # Hyperion used, but it now uses the limits which are computed in C.

    def get_volumes(current_i, current_volume, refined):
        volumes = []
        i = current_i
        for sub in range(8):
            volumes += [current_volume]
            if refined[i]:
                sub_volumes, i = get_volumes(i+1, current_volume / 8., refined)
                volumes += sub_volumes
            else:
                i += 1
        return volumes, i

    volumes_ref, last_i = get_volumes(1, 6 * 5 * 4, o.refined)
    volumes_ref.insert(0, 6 * 5 * 4 * 8.)

    assert_allclose(volumes, volumes_ref)
