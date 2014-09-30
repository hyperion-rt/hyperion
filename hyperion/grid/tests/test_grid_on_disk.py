import h5py

import numpy as np

from ...model.tests.test_helpers import random_id
from .. import CartesianGrid, GridOnDisk


def test_grid_on_disk(tmpdir):

    grid_file = tmpdir.join(random_id()).strpath

    f = h5py.File(grid_file, 'w')

    grid = CartesianGrid([-1., 1.], [-1., 1.], [-1., 1.])
    grid['density'] = []
    grid['density'].append(np.array([[[1.]]]))
    grid.write(f)

    f.close()

    g = GridOnDisk(grid_file)
    assert g['density'].n_dust == 1

