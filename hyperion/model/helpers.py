import numpy as np

from .model import Model
from ..grid import SphericalPolarGrid
from ..util.interpolate import interp1d_fast


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

    # Finalize the model so that the density grids get set
    m = model.to_model()

    # Initialize cumulative optical depth array
    tau_all = np.zeros(m.grid.shape)

    # Loop over envelopes and add cumulative column density
    for i, item in enumerate(m.grid['density']):

        # Find density
        rho = item.array

        # Find optical depth in all cells in radial direction
        dtau = m.grid.widths[0, :, :, :] * rho * m.dust[i].optical_properties.interp_chi_wav(wav)

        # Find cumulative sum starting from the ouside
        tau_all += np.cumsum(dtau[:, :, ::-1], axis=2)

    r = np.zeros(m.grid.shape[:2])
    for ip in range(m.grid.shape[0]):
        for it in range(m.grid.shape[1]):
            tau_col = np.hstack([0., tau_all[ip, it, :]])
            if tau < np.max(tau_col):
                r[ip, it] = interp1d_fast(tau_col, m.grid.r_wall[::-1], tau)
                # print(tau_col, r[ip, it])
            else:
                r[ip, it] = 0.

    return r
