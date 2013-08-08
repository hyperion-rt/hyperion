import numpy as np
import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

# Read in the model
m = ModelOutput('quantity_cartesian.rtout')

# Extract the quantities
g = m.get_quantities()

# Get the wall positions in pc
xw, yw = g.x_wall / pc, g.y_wall / pc

# Make a 2-d grid of the wall positions (used by pcoloarmesh)
X, Y = np.meshgrid(xw, yw)

# Calculate the density-weighted temperature
weighted_temperature = (np.sum(g['temperature'][0].array
                               * g['density'][0].array, axis=2)
                        / np.sum(g['density'][0].array, axis=2))

# Make the plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
c = ax.pcolormesh(X, Y, weighted_temperature)
ax.set_xlim(xw[0], xw[-1])
ax.set_xlim(yw[0], yw[-1])
ax.set_xlabel('x (pc)')
ax.set_ylabel('y (pc)')
cb = fig.colorbar(c)
cb.set_label('Temperature (K)')
fig.savefig('weighted_temperature_cartesian.png', bbox_inches='tight')

# show image

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
c = ax.pcolormesh(X, Y, g['temperature'][0].array[:, 49, :])
ax.set_xlim(xw[0], xw[-1])
ax.set_xlim(yw[0], yw[-1])
ax.set_xlabel('x (pc)')
ax.set_ylabel('y (pc)')
cb = fig.colorbar(c)
cb.set_label('Temperature (K)')
fig.savefig('sliced_temperature_cartesian.png', bbox_inches='tight')
