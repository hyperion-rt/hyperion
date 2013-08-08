
from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import au, lsun, rsun, tsun, msun

# Initialize model
m = AnalyticalYSOModel()

# Set up star
m.star.radius = 1.5 * rsun
m.star.temperature = tsun
m.star.luminosity = lsun

# Set up disk
d = m.add_flared_disk()
d.rmin = 10 * rsun
d.rmax = 30. * au
d.mass = 0.01 * msun
d.p = -1
d.beta = 1.25
d.r_0 = 10. * au
d.h_0 = 0.4 * au
d.dust = 'kmh_lite.hdf5'

# Set up grid
m.set_spherical_polar_grid_auto(400, 100, 1)

# Don't compute temperatures
m.set_n_initial_iterations(0)

# Don't re-emit photons
m.set_kill_on_absorb(True)

# Use raytracing (only important for source here, since no dust emission)
m.set_raytracing(True)

# Compute images using monochromatic radiative transfer
m.set_monochromatic(True, wavelengths=[1.])

# Set up image
i = m.add_peeled_images()
i.set_image_limits(-13 * rsun, 13 * rsun, -13. * rsun, 13 * rsun)
i.set_image_size(256, 256)
i.set_viewing_angles([60.], [20.])
i.set_wavelength_range(1, 1, 1)

# Set number of photons
m.set_n_photons(imaging_sources=10000000, imaging_dust=0,
                raytracing_sources=100000, raytracing_dust=0)

# Write out the model and run it in parallel
m.write('pure_scattering.rtin')
m.run('pure_scattering.rtout', mpi=True)
