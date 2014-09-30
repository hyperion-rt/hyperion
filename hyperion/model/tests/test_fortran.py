from astropy.tests.helper import pytest
import numpy as np

from ...dust import IsotropicDust
from ...util.functions import B_nu

from .. import Model

from .test_helpers import random_id, get_test_dust


def test_point_source_outside_grid(tmpdir):

    dust = get_test_dust()

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), dust)
    m.set_n_photons(initial=100, imaging=0)
    s = m.add_point_source()
    s.position = (-1.5, 0., 0.)
    s.temperature = 5000.
    s.luminosity = 1.
    m.write(tmpdir.join(random_id()).strpath)
    log_file = tmpdir.join(random_id()).strpath
    with pytest.raises(SystemExit) as exc:
        m.run(tmpdir.join(random_id()).strpath, logfile=log_file)
    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert 'photon was not emitted inside a cell' in open(log_file).read()


def test_unsorted_spectrum(tmpdir):

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.set_n_photons(initial=100, imaging=0)
    s = m.add_point_source()
    s._spectrum = {'nu': [3.e20, 2.e10, 1], 'fnu': [1, 2, 3]}
    s.luminosity = 1.
    m.write(tmpdir.join(random_id()).strpath)
    log_file = tmpdir.join(random_id()).strpath
    with pytest.raises(SystemExit) as exc:
        m.run(tmpdir.join(random_id()).strpath, logfile=log_file)
    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert 'spectrum frequency should be monotonically increasing' in open(log_file).read()


def test_spectrum_dust_nooverlap(tmpdir):

    # Set up dust with a narrow frequency range
    nu = np.logspace(8., 10., 100)
    albedo = np.repeat(0.5, 100)
    chi = np.ones(100)
    d = IsotropicDust(nu, albedo, chi)
    d.set_lte_emissivities(10, 0.1, 1000.)

    # Set up model with a source with a wider frequency range
    m = Model()

    s = m.add_point_source()
    s.luminosity = 1.
    nu = np.logspace(5., 12., 1000)
    s.spectrum = (nu, B_nu(nu, 6000.))

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1])

    m.add_density_grid(np.array([[[1.]]]), d)

    m.set_n_photons(initial=1000, imaging=0)

    m.write(tmpdir.join(random_id()).strpath)
    log_file = tmpdir.join(random_id()).strpath
    with pytest.raises(SystemExit) as exc:
        m.run(tmpdir.join(random_id()).strpath, logfile=log_file)
    assert exc.value.args[0] == 'An error occurred, and the run did not ' + \
                                'complete'
    assert 'photon frequency' in open(log_file).read()
    assert 'is outside the range defined' in open(log_file).read()
    assert 'for the dust optical properties' in open(log_file).read()
