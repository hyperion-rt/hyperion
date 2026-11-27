import glob

import numpy as np
import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import au

# Find the output files for all the iterations
filenames = sorted(glob.glob('hseq_?????.rtout'))

fig = plt.figure(figsize=(10, 4))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

previous = None
change = []

for iteration, filename in enumerate(filenames):

    # Extract the grid and the physical quantities
    g = ModelOutput(filename).get_quantities()
    density = g['density'][0].array[0, :, :]
    temperature = g['temperature'][0].array[0, :, :]

    # Find the radial cell closest to 50au
    iw = np.argmin(np.abs(g.w - 50. * au))

    # Plot the vertical density and temperature profiles at that radius
    if iteration in (0, 1, len(filenames) - 1):
        label = 'Iteration {0}'.format(iteration)
        ax1.plot(g.z / au, density[:, iw], label=label)
        ax2.plot(g.z / au, temperature[:, iw], label=label)

    # Find the fraction of the disk mass that has moved compared to the
    # previous iteration
    if previous is not None:
        volumes = g.volumes[0, :, :]
        change.append(np.sum(np.abs(density - previous) * volumes) /
                      np.sum(density * volumes))

    previous = density

for i, c in enumerate(change):
    print('Fraction of mass moved in iteration {0}: {1:.3f}'.format(i + 1, c))

ax1.set_xlim(-50., 50.)
ax1.set_yscale('log')
ax1.set_ylim(1.e-20, 1.e-13)
ax1.set_xlabel('z (au)')
ax1.set_ylabel('Density (g/cm^3)')
ax1.legend(loc='lower center', fontsize='small')

ax2.set_xlim(-50., 50.)
ax2.set_xlabel('z (au)')
ax2.set_ylabel('Temperature (K)')

fig.suptitle('Vertical profiles at w = {0:.0f} au'.format(g.w[iw] / au))
fig.savefig('hseq_profiles.png', bbox_inches='tight')
