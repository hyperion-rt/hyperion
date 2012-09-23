import numpy as np

from .. import Model
from ...util.functions import random_filename
from .test_helpers import get_test_dust


def test_no_initial():
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), get_test_dust())
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=1)
    m.write(random_filename())
    mo = m.run(random_filename())
    g = mo.get_quantities()
    assert 'density' in g