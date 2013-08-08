from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import lsun, rsun, tsun, msun, au

# Initialize model and set up density grid
m = AnalyticalYSOModel()

# Set up the central source
m.star.radius = rsun
m.star.temperature = tsun
m.star.luminosity = lsun

# Set up a simple flared disk
d = m.add_flared_disk()
d.mass = 0.001 * msun
d.rmin = 0.1 * au
d.rmax = 100. * au
d.p = -1
d.beta = 1.25
d.h_0 = 0.01 * au
d.r_0 = au
d.dust = 'kmh_lite.hdf5'

# Specify that the specific energy and density are needed
m.conf.output.output_specific_energy = 'last'
m.conf.output.output_density = 'last'

# Set the number of photons
m.set_n_photons(initial=1000000, imaging=0)

# Set up the grid
m.set_spherical_polar_grid_auto(400, 300, 1)

# Use MRW and PDA
m.set_mrw(True)
m.set_pda(True)

# Write output and run model
m.write('quantity_spherical.rtin')
m.run('quantity_spherical.rtout', mpi=True)
