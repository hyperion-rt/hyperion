import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from scipy.spatial import cKDTree

from hyperion.model import ModelOutput
from hyperion.util.constants import au, pc

# Radius of the cloud and half-width of the domain, as in the setup script
R = 1.e4 * au
L = 1.6 * R

# Pixel positions for a slice through the grid at z=0 (in cm)
n_pix = 350
xs = np.linspace(-L, L, n_pix)
X, Y = np.meshgrid(xs, xs)
pixels = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])

fig, axes = plt.subplots(2, 2, figsize=(9, 9.5))

for column, (name, title) in enumerate([('voronoi_cloud_nowall', 'No wall'),
                                        ('voronoi_cloud_wall', 'Spherical wall')]):

    mo = ModelOutput(name + '.rtout')

    # Read back the Voronoi sites and the density in each cell
    grid = mo.get_quantities()
    sites = np.column_stack([grid.x, grid.y, grid.z])
    rho = grid['density'][0].array

    # Each pixel belongs to the cell of the closest site
    _, cell = cKDTree(sites).query(pixels)
    cell = cell.reshape(X.shape)
    rho_slice = rho[cell]

    if title == 'Spherical wall':
        # With the spherical wall, each cell is cut by the plane tangent to
        # the sphere at the point closest to its site, so we blank out the
        # pixels that lie beyond that plane
        s = sites[cell]
        s_unit = s / np.linalg.norm(s, axis=-1)[:, :, np.newaxis]
        outside = X * s_unit[:, :, 0] + Y * s_unit[:, :, 1] > R
        rho_slice[outside] = np.nan
        cell[outside] = -1

    # Find the boundaries between cells, to overlay them on the density
    edges = np.zeros(X.shape, dtype=bool)
    edges[1:, :] |= cell[1:, :] != cell[:-1, :]
    edges[:, 1:] |= cell[:, 1:] != cell[:, :-1]
    edges_rgba = np.zeros(X.shape + (4,))
    edges_rgba[edges] = (1., 1., 1., 0.4)

    extent = [-L / au / 1000, L / au / 1000, -L / au / 1000, L / au / 1000]

    # Top row: slice through the cells, colored by density
    ax = axes[0, column]
    ax.set_facecolor('black')
    ax.imshow(rho_slice, origin='lower', extent=extent, cmap='magma',
              norm=LogNorm(vmin=5e-21, vmax=3e-18), interpolation='nearest')
    ax.imshow(edges_rgba, origin='lower', extent=extent, interpolation='antialiased')
    ax.add_patch(Circle((0, 0), R / au / 1000, fill=False, color='#7fd3ff',
                        lw=1, ls='--'))
    ax.set_title(title + ': cells at z=0')

    # Bottom row: three-color scattered light image (red=0.8, green=0.55,
    # blue=0.45 microns), with the same arcsinh stretch for both models
    image = mo.get_image(units='MJy/sr', distance=300 * pc)
    rgb = image.val[0][:, :, ::-1]
    rgb = np.clip(np.arcsinh(rgb / 100.) / np.arcsinh(5000. / 100.), 0., 1.)

    ax = axes[1, column]
    ax.imshow(rgb, origin='lower', extent=extent, interpolation='nearest')
    ax.add_patch(Circle((0, 0), R / au / 1000, fill=False, color='#7fd3ff',
                        lw=1, ls='--'))
    ax.set_title(title + ': image at 0.45/0.55/0.8 $\\mu$m')

for ax in axes.ravel():
    ax.set_xlabel('x (1000 au)')
    ax.set_ylabel('y (1000 au)')

fig.tight_layout()
fig.savefig('voronoi_cloud.png', dpi=80)
