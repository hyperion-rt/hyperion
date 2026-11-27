import glob

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from hyperion.model import ModelOutput
from hyperion.util.constants import au

# Find the output files for all the iterations
filenames = sorted(glob.glob('hseq_?????.rtout'))

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)

# Use the same logarithmic color scale for both panels, and show cells with
# zero density or a density below the color scale in the lowest color
norm = LogNorm(vmin=1.e-20, vmax=1.e-12)
cmap = plt.get_cmap('viridis')
cmap = cmap.with_extremes(bad=cmap(0.), under=cmap(0.))

# Show the initial model and the last iteration
for ax, iteration in zip(axes, (0, len(filenames) - 1)):

    # Extract the grid and the density in the (w, z) plane
    g = ModelOutput(filenames[iteration]).get_quantities()
    density = g['density'][0].array[0, :, :]

    # Plot the density using the cell walls, so that the cells are shown
    # with their actual sizes
    image = ax.pcolormesh(g.w_wall / au, g.z_wall / au, density,
                          norm=norm, cmap=cmap)

    # Show a few density contours to make the changes easier to see
    ax.contour(g.w / au, g.z / au, density, levels=[1.e-18, 1.e-16, 1.e-14],
               colors='white', linewidths=0.75)

    ax.set_xlim(0., 100.)
    ax.set_ylim(-50., 50.)
    ax.set_xlabel('w (au)')
    ax.set_title('Iteration {0}'.format(iteration))

axes[0].set_ylabel('z (au)')
fig.colorbar(image, ax=axes, label='Density (g/cm^3)', extend='both')
fig.savefig('hseq_density_maps.png', bbox_inches='tight')
