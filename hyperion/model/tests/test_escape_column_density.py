import numpy as np
import pytest

from .. import Model
from .test_helpers import random_id
from ...dust import IsotropicDust

# The models in this file are uniform-density spheres of radius 1 divided into
# radial shells, for which the minimum column density to escape from anywhere
# in a shell is reached by going radially outwards from its outer wall, so
# that it is equal to the density times (1 - r_outer).

N_SHELLS = 50
RHO = 5.
CHI = 1.


def sphere_model(peeloff_tau_max=None, escape=True):

    m = Model()
    m.set_spherical_polar_grid(np.linspace(0., 1., N_SHELLS + 1), [0., np.pi], [0., 2. * np.pi])

    dust = IsotropicDust([3.e9, 3.e16], [0., 0.], [CHI, CHI])
    dust.set_lte_emissivities(n_temp=10, temp_min=0.1, temp_max=1.e4)

    m.add_density_grid(np.ones((1, 1, N_SHELLS)) * RHO, dust)

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([90.], [0.])
    i.set_track_origin('basic')

    m.set_n_initial_iterations(0)
    m.set_monochromatic(True, wavelengths=[1.])

    m.set_n_photons(imaging_sources=1000, imaging_dust=0)

    if peeloff_tau_max is not None:
        m.set_max_peeloff_optical_depth(peeloff_tau_max)

    if escape:
        m.set_escape_column_density(True)

    return m


def run(tmpdir, m):
    m.write(tmpdir.join(random_id()).strpath)
    return m.run(tmpdir.join(random_id()).strpath)


@pytest.mark.requires_hyperion_binaries
def test_escape_column_density_sphere(tmpdir):

    out = run(tmpdir, sphere_model())

    column = out.get_quantities()['escape_column_density'][0].array[0, 0, :]

    r_outer = np.linspace(0., 1., N_SHELLS + 1)[1:]
    expected = RHO * (1. - r_outer)

    # The estimate can never be lower than the true minimum, and with the
    # default number of rays should be very close to it
    assert np.all(column >= expected - 1e-10)
    np.testing.assert_allclose(column, expected, rtol=0.01, atol=1e-10)

    # The outermost shell is on the edge of the grid
    assert column[-1] == 0.


@pytest.mark.requires_hyperion_binaries
def test_escape_column_density_not_computed(tmpdir):
    out = run(tmpdir, sphere_model(escape=False))
    assert 'escape_column_density' not in out.get_quantities()


def source_flux(tmpdir, **kwargs):
    out = run(tmpdir, sphere_model(**kwargs))
    return out.get_sed(component='source_emit').val[0, 0, 0]


@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize('escape', [False, True])
def test_peeloff_hidden_cells(tmpdir, escape):

    # The source is in the central cell, from which the minimum optical depth
    # to escape is (1 - 1 / N_SHELLS) * RHO * CHI = 4.9, and the optical depth
    # to the observer is RHO * CHI = 5.

    reference = source_flux(tmpdir, escape=False, peeloff_tau_max=1000.) / np.exp(-5.)

    # The optical depth to escape is lower than the maximum, so the photons
    # are peeled off as usual
    flux = source_flux(tmpdir, escape=escape, peeloff_tau_max=10.)
    np.testing.assert_allclose(flux / reference, np.exp(-5.), rtol=1e-10)

    # The optical depth to escape is higher than the maximum. Without the
    # minimum escape column density, the optical depth integration stops once
    # it exceeds the maximum, and with it, the photons are not peeled off at
    # all.
    flux = source_flux(tmpdir, escape=escape, peeloff_tau_max=2.)
    if escape:
        assert flux == 0.
    else:
        assert np.exp(-2.1) <= flux / reference <= np.exp(-2.) * (1 + 1e-10)

