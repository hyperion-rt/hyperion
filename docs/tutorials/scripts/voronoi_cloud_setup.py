import numpy as np

from hyperion.model import Model
from hyperion.util.constants import au, lsun

# Use a fixed random seed so that the cloud is the same every time
rng = np.random.default_rng(2026)

# Radius of the spherical cloud
R = 1.e4 * au

# Background density (in g/cm^3), giving an optical depth of order unity from
# the center to the edge of the cloud at optical wavelengths, not counting the
# clumps
RHO_0 = 3.e-20

# The cloud contains a number of Gaussian clumps, with random positions
# inside the sphere, and a range of sizes and peak densities
n_clumps = 45
directions = rng.normal(size=(n_clumps, 3))
directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]
radii = 0.85 * R * rng.uniform(0., 1., n_clumps) ** (1. / 3.)
clump_centers = directions * radii[:, np.newaxis]
clump_sizes = R * 10 ** rng.uniform(-1.6, -0.9, n_clumps)
clump_amplitudes = 10 ** rng.uniform(1., 2., n_clumps)


def density(x, y, z):
    """
    Return the density at positions x, y, z (in cm).
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # Background density, decreasing from the center to the edge
    rho = RHO_0 * (1. - 0.7 * r / R)
    # Add the clumps
    for center, size, amplitude in zip(clump_centers, clump_sizes, clump_amplitudes):
        d2 = (x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2
        rho = rho + RHO_0 * amplitude * np.exp(-0.5 * d2 / size ** 2)
    return rho


def random_points_in_sphere(n, radius):
    """
    Return n random positions distributed uniformly inside a sphere.
    """
    directions = rng.normal(size=(n, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]
    radii = radius * rng.uniform(0., 1., n) ** (1. / 3.)
    return directions * radii[:, np.newaxis]


# Choose the positions of the Voronoi sites. If you have particles from a
# simulation (e.g. SPH), you can simply use the particle positions. Here we
# place the sites so that they follow the density: a uniform population of
# sites across the whole cloud, plus many more sites inside the clumps, so
# that the cells are smallest where the density is highest.
sites = [random_points_in_sphere(3000, R)]
for center, size, amplitude in zip(clump_centers, clump_sizes, clump_amplitudes):
    n = int(40 * np.sqrt(amplitude))
    sites.append(center + 1.5 * size * rng.normal(size=(n, 3)))
sites = np.vstack(sites)

# All sites must be inside the sphere when using a spherical wall, so we
# drop the few clump sites that fall outside
sites = sites[np.linalg.norm(sites, axis=1) < 0.995 * R]
x, y, z = sites[:, 0], sites[:, 1], sites[:, 2]
print(f"Using {len(x)} Voronoi sites")

# The density in each cell is taken to be the density at the position of its
# site
rho = density(x, y, z)


def setup_model(wall):
    """
    Set up the model, with or without a spherical wall around the cloud.
    """

    m = Model()

    # Set up the Voronoi grid. The cells fill the rectangular domain given by
    # the limits (in cm). If a spherical wall is used, the cells are in
    # addition cut by a sphere, given by its center and radius.
    if wall:
        m.set_voronoi_grid(x, y, z,
                           xmin=-1.6 * R, xmax=1.6 * R,
                           ymin=-1.6 * R, ymax=1.6 * R,
                           zmin=-1.6 * R, zmax=1.6 * R,
                           wall='sphere', wall_args=(0., 0., 0., R))
    else:
        m.set_voronoi_grid(x, y, z,
                           xmin=-1.6 * R, xmax=1.6 * R,
                           ymin=-1.6 * R, ymax=1.6 * R,
                           zmin=-1.6 * R, zmax=1.6 * R)

    # Add the density, with one value per site, in the same order as the sites
    m.add_density_grid(rho, 'kmh_lite.hdf5')

    # Add a hot star at the center and two cooler stars
    for position, luminosity, temperature in [((0., 0., 0.), 3.e4, 20000.),
                                              ((0.45 * R, -0.25 * R, 0.1 * R), 2.e3, 6000.),
                                              ((-0.35 * R, 0.4 * R, -0.2 * R), 1.e3, 4000.)]:
        source = m.add_point_source()
        source.position = position
        source.luminosity = luminosity * lsun
        source.temperature = temperature

    # Compute an image covering the whole domain, viewed down the z axis
    image = m.add_peeled_images(sed=False, image=True)
    image.set_viewing_angles([0.], [0.])
    image.set_image_size(150, 150)
    image.set_image_limits(-1.6 * R, 1.6 * R, -1.6 * R, 1.6 * R)

    # We only compute scattered light images at three wavelengths, so we use
    # monochromatic radiative transfer and do not compute temperatures. In
    # monochromatic mode, photons lose energy each time they scatter, and
    # we stop following them once their energy has decreased by a factor
    # of 10^4.
    m.set_monochromatic(True, wavelengths=[0.45, 0.55, 0.8], energy_threshold=1.e-4)
    m.set_n_initial_iterations(0)
    m.set_n_photons(imaging_sources=500000, imaging_dust=0)

    return m


for name, wall in [('voronoi_cloud_nowall', False),
                   ('voronoi_cloud_wall', True)]:
    m = setup_model(wall)
    m.write(name + '.rtin')
    m.run(name + '.rtout', overwrite=True)
