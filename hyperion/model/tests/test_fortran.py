import pytest
import numpy as np

from .. import Model
from .test_helpers import random_filename, get_test_dust


def test_point_source_outside_grid():

    dust = get_test_dust()

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), dust)
    m.set_n_photons(initial=100, imaging=0)
    s = m.add_point_source()
    s.position = (-1.5, 0., 0.)
    s.temperature = 5000.
    s.luminosity = 1.
    m.write(random_filename())
    log_file = random_filename()
    with pytest.raises(SystemExit) as exc:
        m.run(random_filename(), logfile=log_file)
    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert 'photon was not emitted inside a cell' in open(log_file).read()


def test_unsorted_spectrum():

    dust = get_test_dust()

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.set_n_photons(initial=100, imaging=0)
    s = m.add_point_source()
    s._spectrum = {'nu': [3.e20,2.e10,1], 'fnu': [1,2,3]}
    s.luminosity = 1.
    m.write(random_filename())
    log_file = random_filename()
    with pytest.raises(SystemExit) as exc:
        m.run(random_filename(), logfile=log_file)
    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert 'spectrum frequency should be monotonically increasing' in open(log_file).read()
