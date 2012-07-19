import numpy as np

from ..model import Model
from .test_helpers import random_filename, get_test_dust

def test_mono_zero_prob():

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

    m.write(random_filename())
    m.run(random_filename())
