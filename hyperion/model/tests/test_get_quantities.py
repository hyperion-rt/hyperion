import pytest

import numpy as np

from .. import Model
from ...util.functions import random_id
from .test_helpers import get_test_dust


@pytest.mark.requires_hyperion_binaries
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


@pytest.mark.requires_hyperion_binaries
def test_get_available_components_zero_based(tmpdir):
    # Regression test for get_available_components using 1-based iteration
    # numbering while get_quantities is zero-based
    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), get_test_dust())
    source = m.add_point_source()
    source.luminosity = 1.
    source.temperature = 6000.
    m.set_n_initial_iterations(1)
    m.set_n_photons(initial=1, imaging=0)
    m.write(tmpdir.join(random_id()).strpath)
    mo = m.run(tmpdir.join(random_id()).strpath)
    components = mo.get_available_components(iteration=0)
    assert 'specific_energy' in components
    assert mo.get_available_components() == components
