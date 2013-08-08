import random
random.seed('hyperion')  # ensure that random numbers are the same every time

import numpy as np
from hyperion.model import Model
from hyperion.util.constants import pc, lsun

# Define cell walls
x = np.linspace(-10., 10., 101) * pc
y = np.linspace(-10., 10., 101) * pc
z = np.linspace(-10., 10., 101) * pc

# Initialize model and set up density grid
m = Model()
m.set_cartesian_grid(x, y, z)
m.add_density_grid(np.ones((100, 100, 100)) * 1.e-20, 'kmh_lite.hdf5')

# Generate random sources
for i in range(100):
    s = m.add_point_source()
    xs = random.uniform(-10., 10.) * pc
    ys = random.uniform(-10., 10.) * pc
    zs = random.uniform(-10., 10.) * pc
    s.position = (xs, ys, zs)
    s.luminosity = 10. ** random.uniform(0., 3.) * lsun
    s.temperature = random.uniform(3000., 8000.)

# Specify that the specific energy and density are needed
m.conf.output.output_specific_energy = 'last'
m.conf.output.output_density = 'last'

# Set the number of photons
m.set_n_photons(initial=10000000, imaging=0)

# Write output and run model
m.write('quantity_cartesian.rtin')
m.run('quantity_cartesian.rtout', mpi=True)
