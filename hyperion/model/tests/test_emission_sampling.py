import numpy as np
import pytest

from .. import Model, AnalyticalYSOModel
from .test_helpers import random_id, get_test_dust
from ...dust import IsotropicDust
from ...util.constants import au, lsun, msun, rsun, tsun
from ...util.functions import virtual_file


def sphere_flux(tmpdir, emission_sampling_exponent):

    # A purely absorbing sphere of uniform density and specific energy divided
    # into radial shells, so that the emission of the cells varies by orders
    # of magnitude (in proportion to their volume) and the photons emitted
    # from the inner shells are attenuated by up to exp(-5). Dust photons are
    # emitted at a single wavelength and peeled off, and the flux should not
    # depend on how the cells are sampled.

    m = Model()
    m.set_spherical_polar_grid(np.linspace(0., 1., 51), [0., np.pi], [0., 2. * np.pi])

    dust = IsotropicDust([3.e9, 3.e16], [0., 0.], [1., 1.])
    dust.set_lte_emissivities(n_temp=10, temp_min=0.1, temp_max=1.e4)
    density = np.ones((1, 1, 50)) * 5.
    m.add_density_grid(density, dust, specific_energy=density * 0. + 1.e10)

    s = m.add_point_source()
    s.luminosity = 1.
    s.temperature = 6000.

    i = m.add_peeled_images(sed=True, image=False)
    i.set_viewing_angles([90.], [0.])
    i.set_uncertainties(True)

    m.set_n_initial_iterations(0)
    if emission_sampling_exponent is None:
        m.set_monochromatic(True, wavelengths=[1.])
    else:
        m.set_monochromatic(True, wavelengths=[1.], emission_sampling_exponent=emission_sampling_exponent)
    m.set_n_photons(imaging_sources=0, imaging_dust=20000)
    m.set_seed(-12345)

    m.write(tmpdir.join(random_id()).strpath)
    out = m.run(tmpdir.join(random_id()).strpath)
    sed = out.get_sed(uncertainties=True)
    return sed.val[0, 0, 0], sed.unc[0, 0, 0]


@pytest.mark.requires_hyperion_binaries
def test_emission_sampling_exponent_default(tmpdir):
    # An exponent of 1 is the default and should give identical results
    assert sphere_flux(tmpdir, 1.) == sphere_flux(tmpdir, None)


@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize('exponent', [0., 0.5])
def test_emission_sampling_exponent_unbiased(tmpdir, exponent):
    # Other exponents give different photons but should agree with the
    # default within the uncertainties
    flux_ref, unc_ref = sphere_flux(tmpdir, 1.)
    flux, unc = sphere_flux(tmpdir, exponent)
    assert flux != flux_ref
    assert abs(flux - flux_ref) < 4 * np.hypot(unc, unc_ref)
    assert unc / flux < 0.05


@pytest.mark.parametrize('value', [0, 0.25, 0.5, 1, np.float32(0.5)])
def test_emission_sampling_exponent_valid(value):
    m = Model()
    m.set_monochromatic(True, wavelengths=[1.], emission_sampling_exponent=value)
    assert m._monochromatic_emission_sampling_exponent == value


@pytest.mark.parametrize(('value', 'exception', 'message'), [
    ('a', TypeError, "emission_sampling_exponent should be a scalar value"),
    ([0.5, 0.5], TypeError, "emission_sampling_exponent should be a scalar value"),
    (-0.1, ValueError, "emission_sampling_exponent should be between 0 and 1"),
    (1.5, ValueError, "emission_sampling_exponent should be between 0 and 1"),
    (float('nan'), ValueError, "emission_sampling_exponent should be between 0 and 1"),
])
def test_emission_sampling_exponent_invalid(value, exception, message):
    m = Model()
    with pytest.raises(exception) as exc:
        m.set_monochromatic(True, wavelengths=[1.], emission_sampling_exponent=value)
    assert exc.value.args[0] == message


def test_emission_sampling_exponent_requires_monochromatic():
    m = Model()
    with pytest.raises(ValueError) as exc:
        m.set_monochromatic(False, emission_sampling_exponent=0.5)
    assert exc.value.args[0] == "emission_sampling_exponent can only be set in monochromatic mode"
    # The default value is allowed
    m.set_monochromatic(False, emission_sampling_exponent=1.)


@pytest.mark.parametrize('value', [1., 0.5, 0.])
def test_emission_sampling_exponent_io(value):
    m1 = Model()
    m1.set_monochromatic(True, wavelengths=[1.], emission_sampling_exponent=value)
    v = virtual_file()
    m1._write_monochromatic(v)
    # The default is not written out, so that input files are unchanged
    assert ('emission_sampling_exponent' in v.attrs) is (value != 1.)
    m2 = Model()
    m2.set_monochromatic(True, wavelengths=[2.], emission_sampling_exponent=0.3)
    m2._read_monochromatic(v)
    assert m2._monochromatic_emission_sampling_exponent == value


def test_emission_sampling_exponent_analytical_yso():
    # The setting should be preserved when converting an AnalyticalYSOModel
    # to a Model, which is done when writing it out
    m = AnalyticalYSOModel()
    m.star.radius = rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun
    disk = m.add_flared_disk()
    disk.mass = 1.e-2 * msun
    disk.rmin = 0.1 * au
    disk.r_0 = au
    disk.rmax = 300. * au
    disk.h_0 = 0.01 * au
    disk.p = -1.
    disk.beta = 1.25
    disk.dust = get_test_dust()
    m.set_spherical_polar_grid_auto(10, 10, 1)
    m.set_n_initial_iterations(0)
    m.set_monochromatic(True, wavelengths=[1.], emission_sampling_exponent=0.5)
    m.set_n_photons(imaging_sources=10, imaging_dust=10)
    assert m.to_model()._monochromatic_emission_sampling_exponent == 0.5
