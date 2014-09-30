import h5py

import numpy as np

from numpy.testing import assert_allclose

from ..model import Model
from .test_helpers import random_id, get_test_dust
from ...grid import CartesianGrid


def test_use_grid_from_file(tmpdir):

    grid_file = tmpdir.join(random_id()).strpath
    input_file = tmpdir.join(random_id()).strpath
    output_file = tmpdir.join(random_id()).strpath

    f = h5py.File(grid_file, 'w')

    grid = CartesianGrid([-1., 1.], [-1., 1.], [-1., 1.])
    grid['density'] = []
    grid['density'].append(np.array([[[1.]]]))
    grid.write(f)

    f.close()

    model = Model()

    model.use_grid_from_file(grid_file, dust=[get_test_dust()])

    model.set_n_photons(initial=1, imaging=0)
    model.set_n_initial_iterations(3)

    source = model.add_point_source()
    source.luminosity = 1.
    source.temperature = 1000.

    model.write(input_file)
    model.run(output_file)

