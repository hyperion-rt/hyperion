import numpy as np

from hyperion.model import Model
from hyperion.util.constants import pc, c

# The following value is taken from Mathis, Mezger, and Panagia (1983)
FOUR_PI_JNU = 0.0217

# Initialize model
m = Model()

# Set up grid
m.set_spherical_polar_grid([0., 1.001 * pc],
                           [0., np.pi],
                           [0., 2. * np.pi])

# Read in MMP83 spectrum
wav, jlambda = np.loadtxt('mmp83.txt', unpack=True)
nu = c / (wav * 1.e-4)
jnu = jlambda * wav / nu

# Set up the source - note that the normalization of the spectrum is not
# important - the luminosity is set separately.
s = m.add_external_spherical_source()
s.radius = pc
s.spectrum = (nu, jnu)
s.luminosity = np.pi * pc * pc * FOUR_PI_JNU

# Add an inside observer with an all-sky camera
sed = m.add_peeled_images(sed=False, image=True)
sed.set_inside_observer((0., 0., 0.))
sed.set_image_limits(180., -180., -90., 90.)
sed.set_image_size(256, 128)
sed.set_wavelength_range(100, 0.01, 1000.)

# Use raytracing for high signal-to-noise
m.set_raytracing(True)

# Don't compute the temperature
m.set_n_initial_iterations(0)

# Only include photons from the source (since there is no dust)
m.set_n_photons(imaging=0,
                raytracing_sources=10000000,
                raytracing_dust=0)

# Write out and run the model
m.write('example_isrf.rtin')
m.run('example_isrf.rtout')
