import glob

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from hyperion.model import ModelOutput
from hyperion.util.constants import au

# Find the output files for all the iterations
filenames = sorted(glob.glob('hseq_?????.rtout'))

# Top row: the whole disk, bottom row: the inner edge of the disk
fig, axes = plt.subplots(2, 2, figsize=(10, 8.5))

# Use the same logarithmic color scale for all panels, and show cells with
# zero density or a density below the color scale in the lowest color
norm = LogNorm(vmin=1.e-20, vmax=1.e-7)
cmap = plt.get_cmap('viridis')
cmap = cmap.with_extremes(bad=cmap(0.), under=cmap(0.))

# Show the initial model and the last iteration
for column, iteration in enumerate((0, len(filenames) - 1)):

    # Extract the grid and the density in the (w, z) plane
    g = ModelOutput(filenames[iteration]).get_quantities()
    density = g['density'][0].array[0, :, :]

    for row in (0, 1):

        ax = axes[row, column]

        # Plot the density using the cell walls, so that the cells are shown
        # with their actual sizes
        image = ax.pcolormesh(g.w_wall / au, g.z_wall / au, density,
                              norm=norm, cmap=cmap)

        # Show a few density contours to make the changes easier to see,
        # ignoring the empty cells inside the inner radius of the disk
        ax.contour(g.w / au, g.z / au, np.ma.masked_equal(density, 0.),
                   levels=[1.e-17, 1.e-14, 1.e-11, 1.e-8],
                   colors='white', linewidths=0.75)

        if row == 0:
            # The radial cells are spaced logarithmically, so show w on a log
            # scale, on which all the cells have the same width (the first
            # wall is at w=0)
            ax.set_xscale('log')
            ax.set_xlim(g.w_wall[1] / au, 10.)
            ax.set_ylim(-6., 6.)
            ax.set_title('Iteration {0}'.format(iteration))
        else:
            # Zoom in on the inner edge of the disk
            ax.set_xlim(0., 0.1)
            ax.set_ylim(-0.03, 0.03)

        ax.set_xlabel('w (au)')
        if column == 0:
            ax.set_ylabel('z (au)')

fig.colorbar(image, ax=axes, label='Density (g/cm^3)', extend='both')
fig.savefig('hseq_density_maps.png', bbox_inches='tight')
