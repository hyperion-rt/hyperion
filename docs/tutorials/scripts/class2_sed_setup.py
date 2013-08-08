import numpy as np

from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import rsun, au, msun, sigma

# Initalize the model
m = AnalyticalYSOModel()

# Set the stellar parameters
m.star.radius = 2. * rsun
m.star.temperature = 4000.
m.star.luminosity = 4 * (2. * rsun) ** 2 * sigma * 4000 ** 4

# Add a flared disk
disk = m.add_flared_disk()
disk.mass = 0.01 * msun
disk.rmin = 10 * m.star.radius
disk.rmax = 200 * au
disk.r_0 = m.star.radius
disk.h_0 = 0.01 * disk.r_0
disk.p = -1.0
disk.beta = 1.25
disk.dust = 'kmh_lite.hdf5'

# Use raytracing to improve s/n of thermal/source emission
m.set_raytracing(True)

# Use the modified random walk
m.set_mrw(True, gamma=2.)

# Set up grid
m.set_spherical_polar_grid_auto(399, 199, 1)

# Set up SED for 10 viewing angles
sed = m.add_peeled_images(sed=True, image=False)
sed.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))
sed.set_wavelength_range(150, 0.02, 2000.)
sed.set_track_origin('basic')

# Set number of photons
m.set_n_photons(initial=1e5, imaging=1e6,
                raytracing_sources=1e4, raytracing_dust=1e6)

# Set number of temperature iterations
m.set_n_initial_iterations(5)

# Write out file
m.write('class2_sed.rtin')
m.run('class2_sed.rtout', mpi=True)
