from hyperion.model import AnalyticalYSOModel
from hyperion.model.helpers import run_with_vertical_hseq
from hyperion.util.constants import au, lsun, rsun, tsun, msun

# Initialize model
m = AnalyticalYSOModel()

# Set up star (the stellar mass is required to compute the hydrostatic
# equilibrium)
m.star.radius = 2. * rsun
m.star.temperature = 4000.
m.star.luminosity = 5 * lsun
m.star.mass = 0.5 * msun

# Set up disk
d = m.add_flared_disk()
d.rmin = 10 * rsun
d.rmax = 100. * au
d.mass = 0.01 * msun
d.p = -1
d.beta = 1.25
d.r_0 = 100. * au
d.h_0 = 10. * au
d.dust = 'kmh_lite.hdf5'

# Set up a cylindrical polar grid
m.set_cylindrical_polar_grid_auto(100, 100, 1)

# Set up the temperature calculation, using the modified random walk and the
# partial diffusion approximation to speed up the calculation in the optically
# thick midplane
m.set_n_initial_iterations(5)
m.set_n_photons(initial=1000000, imaging=0)
m.set_mrw(True)
m.set_pda(True)

# Iterate the radiative transfer and the vertical density structure 10 times
run_with_vertical_hseq('hseq', m, n_iter=10, overwrite=True)
