import numpy as np

from ..model import Model
from .test_helpers import random_id, get_test_dust


def test_mono_zero_prob(tmpdir):

    # Check that when total probability is zero in a given dust type for a given wavelength, the code doesn't crash

    dust = get_test_dust()

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.]]]), dust)
    m.add_density_grid(np.array([[[0.5]]]), dust, merge_if_possible=False)

    image = m.add_peeled_images(sed=True, image=True)
    image.set_image_limits(-2., 2., -2., 2.)
    image.set_image_size(20, 20)
    image.set_viewing_angles([45.], [45.])

    m.set_minimum_temperature(10.)

    m.set_monochromatic(True, wavelengths=[0.01, 0.1, 1., 10., 100., 1000.])

    m.set_n_initial_iterations(0)

    m.set_n_photons(imaging_sources=0, imaging_dust=100)

    m.write(tmpdir.join(random_id()).strpath)
    m.run(tmpdir.join(random_id()).strpath)


def test_check_weighting(tmpdir):
    '''
    This is a regression test for a bug that caused incorrect weighting of the
    monochromatic fluxes in the presence of multiple dust populations with one
    dust population having zero mean_prob at a given frequency
    '''

    d = get_test_dust()

    # First model, two dust types, very different energies, so part of the SED
    # comes from only one source

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.e-10]]]), d, specific_energy=np.array([[[1.e8]]]))
    m.add_density_grid(np.array([[[1.e-10]]]), d, specific_energy=np.array([[[1.e-4]]]))

    image = m.add_peeled_images(sed=True, image=False)
    image.set_viewing_angles([45.], [45.])
    image.set_track_origin('detailed')

    m.set_monochromatic(True, wavelengths=np.logspace(-1., 4., 10))

    m.set_n_initial_iterations(0)

    m.set_n_photons(imaging_sources=0, imaging_dust=10000)

    m.write(tmpdir.join(random_id()).strpath)
    mo1 = m.run(tmpdir.join(random_id()).strpath)

    # Same model but without second dust population. Since the model is
    # optically thin, the SEDs for the first dust population should agree.

    m = Model()
    m.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    m.add_density_grid(np.array([[[1.e-10]]]), d, specific_energy=np.array([[[1.e8]]]))

    image = m.add_peeled_images(sed=True, image=False)
    image.set_viewing_angles([45.], [45.])
    image.set_track_origin('detailed')

    m.set_monochromatic(True, wavelengths=np.logspace(-1., 4., 10))

    m.set_n_initial_iterations(0)

    m.set_n_photons(imaging_sources=0, imaging_dust=10000)

    m.write(tmpdir.join(random_id()).strpath)
    mo2 = m.run(tmpdir.join(random_id()).strpath)

    _, nufnu1 = mo1.get_sed(inclination=-1, aperture=-1, component='dust_emit', dust_id=0)
    _, nufnu2 = mo2.get_sed(inclination=-1, aperture=-1, component='dust_emit', dust_id=0)

    THRESHOLD = 1.02

    # Need to ignore first element since it is zero
    assert np.all((nufnu1[1:] / nufnu2[1:] < THRESHOLD) & (nufnu2[1:] / nufnu1[1:] < THRESHOLD))
