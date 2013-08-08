import numpy as np

from hyperion.model import Model
from hyperion.util.constants import pc, lsun

# Initialize model
m = Model()

# Set up 64x64x64 cartesian grid
w = np.linspace(-pc, pc, 64)
m.set_cartesian_grid(w, w, w)

# Add density grid with constant density and add a higher density cube inside to
# cause a shadow.
density = np.ones(m.grid.shape) * 1e-21
density[26:38, 26:38, 26:38] = 1.e-18
m.add_density_grid(density, 'kmh_lite.hdf5')

# Add a point source in the center
s = m.add_point_source()
s.position = (0.4 * pc, 0., 0.)
s.luminosity = 1000 * lsun
s.temperature = 6000.

# Add multi-wavelength image for a single viewing angle
image = m.add_peeled_images(sed=False, image=True)
image.set_wavelength_range(20, 1., 1000.)
image.set_viewing_angles([60.], [80.])
image.set_image_size(400, 400)
image.set_image_limits(-1.5 * pc, 1.5 * pc, -1.5 * pc, 1.5 * pc)

# Set runtime parameters
m.set_n_initial_iterations(5)
m.set_raytracing(True)
m.set_n_photons(initial=4e6, imaging=4e7,
                raytracing_sources=1, raytracing_dust=1e7)

# Write out input file
m.write('simple_cube.rtin')
m.run('simple_cube.rtout', mpi=True)
