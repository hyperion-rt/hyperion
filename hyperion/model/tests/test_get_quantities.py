import numpy as np

from .. import Model
from ...util.functions import random_id
from .test_helpers import get_test_dust


def test_no_initial(tmpdir):
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), get_test_dust())
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging=1)
    m.write(tmpdir.join(random_id()).strpath)
    mo = m.run(tmpdir.join(random_id()).strpath)
    g = mo.get_quantities()
    assert 'density' in g
