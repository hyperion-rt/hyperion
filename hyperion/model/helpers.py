import multiprocessing
import numpy as np

from .model import Model
from ..grid import SphericalPolarGrid
from ..util.interpolate import interp1d_fast
from .analytical_yso_model import AnalyticalYSOModel


def find_last_iteration(file_handle):
    max_iteration = 0
    for group_name in file_handle:
        if "iteration" in group_name:
            iteration = int(group_name.split('_')[1])
            max_iteration = max(iteration, max_iteration)
    return max_iteration


def tau_to_radius(model, tau, wav):
    """
    Given a Model instance with a spherical polar coordinate grid, find the
    radius from which the optical depth to escape radially is a fixed value.

    This only works for spherical polar grids, but works for 1-, 2-, and
    3-d grids.

    Parameters
    ----------
    model : `~hyperion.model.Model` instance
    tau : float
        The optical depth for which to find the surface
    wav : float
        The wavelength at which the optical depth is defined

    Returns
    -------
    r : np.ndarray
        The radius or radii at which the optical depth to escape radially
        is ``tau`` at ``wav``. This is a scalar, a 1-d, or a 2-d array
        depending on the dimensionality of the grid.
    """

    if not isinstance(model, Model):
        raise TypeError("model should be a Model instance")

    if model.grid is None:
        raise Exception("Coordinate grid has not been defined")

    if not isinstance(model.grid, SphericalPolarGrid):
        raise TypeError("This method can only be called for spherical polar grids")

    # Initialize cumulative optical depth array
    tau_all = np.zeros(model.grid.shape)

    # Loop over envelopes and add cumulative column density
    for i, item in enumerate(model.grid['density']):

        # Find density
        rho = item.array

        # Find optical depth in all cells in radial direction
        dtau = model.grid.widths[0, :, :, :] * rho * model.dust[i].optical_properties.interp_chi_wav(wav)

        # Find cumulative sum starting from the ouside
        tau_all += np.cumsum(dtau[:, :, ::-1], axis=2)

    r = np.zeros(model.grid.shape[:2])
    for ip in range(model.grid.shape[0]):
        for it in range(model.grid.shape[1]):
            tau_col = np.hstack([0., tau_all[ip, it, :]])
            if tau < np.max(tau_col):
                r[ip, it] = interp1d_fast(tau_col, model.grid.r_wall[::-1], tau)
            else:
                r[ip, it] = 0.

    return r


def hseq_profile(w, z, temperature, mstar, mu=2.279):
    """
    Compute the new (normalized) density profile
    corresponding to a given temperature profile

    Parameters
    ----------
    w : float
        The cylindrical radius at which to compute the profile (in cm)
    z : np.ndarray
        The z coordinates of the cells at radius w (in cm)
    temperature : np.ndarray
        The temperatures in the cells (in K)
    mstar : float
        The mass of the star (in g)
    """

    from hyperion.util.constants import G, m_h, k
    from ..util.integrate import integrate, integrate_subset

    # Compute the integrand
    integrand = z / temperature / (w ** 2 + z ** 2) ** 1.5

    # Compute the integral for all cells
    # TODO - inefficient to compute integral from scratch - optimize
    i = np.array([integrate_subset(z, integrand, 0., zmax) for zmax in z])
    i[z < 0] = -i[z < 0]

    # Compute the factor for the integrand
    factor = G * mstar * mu * m_h / k

    # Compute the profile
    density = np.exp(-i * factor) / temperature

    # Normalize the density profile
    density = density / integrate(z, density)

    return density


# The mean molecular weight of H2 + He is given by:
#
# mu = 4 * (X + 1) / (X + 2)
#
# where X is the mass fraction of Helium to Hydrogren. Assuming
#
# X = 0.32476319350473615
#
# gives:
#
# mu = 2.279


def run_with_vertical_hseq(prefix, model, n_iter=10, mpi=False,
                           n_processes=multiprocessing.cpu_count(),
                           overwrite=False):
    """
    Run a model with vertical hydrostatic equilibrium.

    .. note:: this is an experimental function that is currently in
              development. Please use with care!

    The hydrostatic equilibrium condition is only applied to the disk
    components. The following requirements must be met:

    - The model should be an AnalyticalYSOModel
    - The model should be defined on a cylindrical polar grid
    - The stellar mass should be set
    - The model should include at least one disk

    The dust properties for the model can be specified as dust or dust+gas
    densities as this does not have an impact on this calculation - however,
    the hydrostatic equilibrium is computed assuming an H2 + He mix of gas
    (i.e. mu=2.279). Note that this calculation also ignores the effects of
    self-gravity in the disk, which might be important for more massive disks.

    Parameters
    ----------
    prefix : str
        The prefix for the output
    model : `~hyperion.model.analytical_yso_model.AnalyticalYSOModel`
        The model to run
    n_iter : int, optional
        The number of iterations to run the model for
    mpi : bool, optional
        Whether to run the model in parallel
    n_processes : int, optional
        The number of processes to use if ``mpi`` is ``True``
    overwrite : bool, optional
        Whether to overwrite previous files
    """

    from ..grid import CylindricalPolarGrid
    from .model_output import ModelOutput
    from ..util.integrate import integrate

    if not isinstance(model, AnalyticalYSOModel):
        raise TypeError("Can only run hydrostatic equilibrium for AnalyticalYSOModel instances")

    if model.grid['grid_type'] != 'cylindrical':
        raise TypeError("Can only run hydrostatic equilibrium for models with cylindrical polar grids")

    if model.star.mass is None:
        raise ValueError("Stellar mass needs to be defined for calculation of hydrostatic equilibrium")

    if len(model.disks) == 0:
        raise ValueError("Can only run hydrostatic equilibrium for models with disks")
    else:
        n_disks = len(model.disks)

    # Write out initial model
    model.write(prefix + '_00000.rtin', overwrite=overwrite, merge_if_possible=False)

    # Run the initial model
    mo = model.run(prefix + '_00000.rtout', overwrite=overwrite,
                   mpi=mpi, n_processes=n_processes)

    previous = prefix + '_00000.rtout'

    for iteration in range(1, n_iter + 1):

        # Read in output
        mo = ModelOutput(previous)

        # Extract the quantities
        g = mo.get_quantities()

        # Get wall positions
        rw, zw = g.w_wall, g.z_wall

        # Make a 2-d grid of wall positions
        R, Z = np.meshgrid(rw, zw)

        # Extract density and temperature
        density = g['density']
        temperature = g['temperature']

        # TODO: need to find a better way than just assuming the first n
        # density grids are disks

        for idisk in range(n_disks):

            # Vertically extrapolate temperatures
            for i in range(len(g.w)):
                for j in range(len(g.p)):
                    reset = temperature[idisk].array[j, :, i] < 1.
                    temperature[idisk].array[j, reset, i] = np.max(temperature[idisk].array[j, :, i])  # shouldn't be max, but will do for now

            # Compute new density
            for i in range(len(g.w)):
                for j in range(len(g.p)):
                    density[idisk].array[j, :, i] = hseq_profile(g.w[i], g.z, temperature[idisk].array[j, :, i], model.star.mass) * integrate(g.z, density[idisk].array[j, :, i])

        # Instantiate new model based on previous
        m = Model.read(previous)

        # Override the density
        m.grid['density'] = density

        # Write and run
        m.write('{0:s}_{1:05d}.rtin'.format(prefix, iteration), overwrite=overwrite)
        m.run('{0:s}_{1:05d}.rtout'.format(prefix, iteration),
              overwrite=overwrite, mpi=mpi, n_processes=n_processes)

        previous = '{0:s}_{1:05d}.rtout'.format(prefix, iteration)
