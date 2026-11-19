import pytest
import numpy as np

from ..model import Model
from .test_helpers import get_realistic_test_dust, random_id
from ...util.constants import au, lsun, rsun, pi

# The following tests ensure that the correct temperatures are returned
# for the MRW, for single and for multiple dust populations.

D_REF = np.logspace(-5., 12., 18)

T_REF = [24.75280,
         24.66414,
         24.52175,
         21.97109,
         15.53059,
         10.76363,
         7.810127,
         6.672520,
         6.902798,
         16.64318,
         53.08394,
         162.0158,
         438.9026,
         1013.141,
         2156.520,
         4642.825,
         9948.065,
         21211.08]


@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize(('density_ref', 'temperature_ref'), list(zip(D_REF, T_REF)))
def test_single_temperature(tmpdir, density_ref, temperature_ref):

    dust = get_realistic_test_dust()

    m = Model()

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    source = m.add_point_source()
    source.luminosity = 1.
    source.temperature = 6000.

    m.add_density_grid(np.ones(m.grid.shape) * density_ref, dust)

    m.set_n_photons(initial=1000, imaging=0)

    m.set_n_initial_iterations(30)
    m.set_convergence(True, percentile=99., absolute=2., relative=1.02)

    m.set_mrw(True, gamma=2)

    m.set_max_interactions(1000000000)

    m.write(tmpdir.join(random_id()).strpath)
    mo = m.run(tmpdir.join(random_id()).strpath)

    grid = mo.get_quantities()
    temperature = grid['temperature'][0].array[0, 0, 0]

    assert temperature_ref / temperature < 1.1 and temperature / temperature_ref < 1.1


@pytest.mark.requires_hyperion_binaries
@pytest.mark.parametrize(('density_ref', 'temperature_ref'), list(zip(D_REF, T_REF)))
def test_multi_temperature(tmpdir, density_ref, temperature_ref):

    dust = get_realistic_test_dust()

    m = Model()

    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])

    source = m.add_point_source()
    source.luminosity = 1.
    source.temperature = 6000.

    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.1, dust)
    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.2, dust)
    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.3, dust)
    m.add_density_grid(np.ones(m.grid.shape) * density_ref * 0.4, dust)

    m.set_n_photons(initial=1000, imaging=0)

    m.set_n_initial_iterations(30)
    m.set_convergence(True, percentile=99., absolute=2., relative=1.02)

    m.set_mrw(True, gamma=2)

    m.set_max_interactions(1000000000)

    m.write(tmpdir.join(random_id()).strpath)
    mo = m.run(tmpdir.join(random_id()).strpath)

    grid = mo.get_quantities()

    for i in range(4):
        temperature = grid['temperature'][i].array[0, 0, 0]
        assert temperature_ref / temperature < 1.1 and temperature / temperature_ref < 1.1


# The following tests check the state of photons that leave the MRW diffusion
# sphere, using the SED of an optically thick spherical shell computed purely
# from peeled-off photons (no raytracing).


@pytest.fixture(scope='module')
def mrw_shell_seds(tmp_path_factory):

    tmpdir = tmp_path_factory.mktemp('mrw_shell')

    dust = get_realistic_test_dust()

    m = Model()

    # Empty central cavity followed by ten logarithmically spaced shells
    r = np.hstack([0., np.logspace(np.log10(5 * au), np.log10(2000 * au), 11)])
    m.set_spherical_polar_grid(r, [0., pi], [0., 2 * pi])

    # Density falling off as r^-2, normalized to a radial optical depth of
    # 1000 at 0.55 microns, so that the MRW is used in the inner shells
    rho = (0.5 * (r[1:] + r[:-1])) ** -2.
    rho[0] = 0.
    rho *= 1000. / (np.sum(rho * np.diff(r)) * dust.optical_properties.interp_chi_wav(0.55))
    m.add_density_grid(rho[np.newaxis, np.newaxis, :], dust)

    source = m.add_spherical_source()
    source.luminosity = 100 * lsun
    source.radius = 3 * rsun
    source.temperature = 6000.

    m.set_n_initial_iterations(3)
    m.set_n_photons(initial=10000, imaging=10000)
    m.set_raytracing(False)
    m.set_mrw(True, gamma=2.)
    m.set_max_interactions(100000000)
    m.set_seed(-12345)

    sed = m.add_peeled_images(sed=True, image=False)
    sed.set_viewing_angles([45.], [0.])
    sed.set_wavelength_range(40, 0.3, 3000.)
    sed.set_track_origin('basic')

    m.write(str(tmpdir / 'mrw_shell.rtin'))
    mo = m.run(str(tmpdir / 'mrw_shell.rtout'))

    seds = {}
    for component in ['total', 'dust_emit', 'dust_scat']:
        seds[component] = mo.get_sed(inclination=0, aperture=-1, units='ergs/s',
                                     component=component)

    return seds


@pytest.mark.requires_hyperion_binaries
def test_mrw_exit_opacity(mrw_shell_seds):

    # When a photon leaves the MRW diffusion sphere with a new frequency, it
    # needs to be propagated with the opacity for that new frequency. If it
    # keeps the opacity of the frequency it had before the MRW, near-infrared
    # photons leak out of the shell with far-infrared opacities, and the
    # fraction of the flux emerging between 1 and 3 microns, which should be
    # vanishingly small through an optical depth this high, is around 1e-4.

    total = mrw_shell_seds['total']
    near_infrared = (total.wav > 1.) & (total.wav < 3.)

    assert total.val[near_infrared].sum() < 1.e-10 * total.val.sum()


@pytest.mark.requires_hyperion_binaries
def test_mrw_exit_origin(mrw_shell_seds):

    # A photon leaving the MRW diffusion sphere is thermal dust emission, so
    # it should not be flagged as scattered even if its last interaction
    # before the MRW was a scattering. The albedo of the test dust is zero
    # beyond 150 microns, so there can be no genuinely scattered light at
    # these wavelengths.

    dust_emit = mrw_shell_seds['dust_emit']
    dust_scat = mrw_shell_seds['dust_scat']
    far_infrared = dust_emit.wav > 200.

    assert dust_emit.val[far_infrared].sum() > 0.
    assert dust_scat.val[far_infrared].sum() < 1.e-6 * dust_emit.val[far_infrared].sum()
