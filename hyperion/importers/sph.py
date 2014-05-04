import numpy as np


def DEFAULT_STOPPING_CRITERION(x, y, z, dx, dy, dz, px, py, pz, sigma):
    return len(px) <= 2


def refine(x, y, z, dx, dy, dz, px, py, pz, sigma, mass, levels_remaining, stopping_criterion):

    if stopping_criterion(x, y, z, dx, dy, dz, px, py, pz, sigma) or levels_remaining == 0:
        return [False]

    b_all = [True]

    px_pos = px > x
    py_pos = py > y
    pz_pos = pz > z

    for zcomp, zsub in ((~pz_pos, z - dz * 0.5),(pz_pos, z + dz  *0.5)):
        for ycomp, ysub in ((~py_pos, y - dy * 0.5),(py_pos, y + dy  *0.5)):
            for xcomp, xsub in ((~px_pos, x - dx * 0.5),(px_pos, x + dx  *0.5)):

                keep = xcomp & ycomp & zcomp

                b = refine(xsub, ysub, zsub,
                           dx * 0.5, dy * 0.5, dz * 0.5,
                           px[keep], py[keep], pz[keep], sigma[keep], mass[keep],
                           levels_remaining - 1,
                           stopping_criterion)

                b_all += b

    return b_all


def discretize_wrapper(args):
    from ._discretize_sph import _discretize_sph_func
    return _discretize_sph_func(*args)


def construct_octree(x, y, z, dx, dy, dz, px, py, pz, sigma, mass,
                     n_levels=None,
                     stopping_criterion=DEFAULT_STOPPING_CRITERION,
                     mode='exact'):
    """
    Construct an Octree grid from SPH particles.

    This will set up the octree, and then compute the densities. It is
    equivalent to calling
    :func:`~hyperion.importers.sph.compute_octree_geometry` and
    func:`~hyperion.importers.sph.compute_octree_density`.

    Parameters
    ----------
    x, y, z : float
        The center of the top-level cell
    dx, dy, dz : float
        The half-width of the top-level cell
    px, py, pz : np.ndarray
        The positions of the SPH particles
    sigma : np.ndarray
        The sigma of the SPH particles (assuming a Gaussian Kernel)
    mass : np.ndarray
        The mass of the SPH particles
    n_levels : int, optional
        Maximum number of levels in the octree. If not specified, then the
        octree will keep dividing until all cells contain at most one particle.
    stopping_criterion : func, optional
        A function that is used to determine whether to stop refining a cell.
        If not set, then refinement stops once there are two or fewer particles
        in a cell.
    mode : { 'fast', 'exact' }
        Whether to properly compute the integral over the kernel ('exact') or
        simply count the particles in each cell ('fast')

    Returns
    -------
    grid : `~hyperion.grid.octree_grid.OctreeGrid`
        The octree grid with the densities set
    """

    # Set up geometrical OctreeGrid
    octree = compute_octree_geometry(x, y, z, dx, dy, dz, px, py, pz, sigma, mass,
                                     n_levels=n_levels, stopping_criterion=stopping_criterion)

    # Compute densities and add to octree
    density = compute_octree_densities(octree, px, py, pz, sigma, mass, mode=mode)
    octree['density'] = []
    octree['density'].append(density)

    return octree


def compute_octree_geometry(x, y, z, dx, dy, dz, px, py, pz, sigma, mass,
                            n_levels=None,
                            stopping_criterion=DEFAULT_STOPPING_CRITERION):
    """
    Generate an Octree grid from SPH particles.

    Note that this does not set the densities.

    Parameters
    ----------
    x, y, z : float
        The center of the top-level cell
    dx, dy, dz : float
        The half-width of the top-level cell
    px, py, pz : np.ndarray
        The positions of the SPH particles
    sigma : np.ndarray
        The sigma of the SPH particles (assuming a Gaussian Kernel)
    mass : np.ndarray
        The mass of the SPH particles
    n_levels : int, optional
        Maximum number of levels in the octree. If not specified, then the
        octree will keep dividing until all cells contain at most one particle.
    stopping_criterion : func, optional
        A function that is used to determine whether to stop refining a cell.
        If not set, then refinement stops once there are two or fewer particles
        in a cell.

    Returns
    -------
    grid : `~hyperion.grid.octree_grid.OctreeGrid`
        The octree grid
    """

    if n_levels is None:
        n_levels = np.inf

    # Compute refined array
    refined = refine(x, y, z, dx, dy, dz, px, py, pz, sigma,
                     mass, n_levels, stopping_criterion)

    # Set up OctreeGrid instance
    from ..grid import OctreeGrid
    octree = OctreeGrid(x, y, z, dx, dy, dz, refined)

    return octree


def compute_octree_densities(octree, px, py, pz, sigma, mass, mode='exact'):
    """
    Compute the density of the cells in an Octree.

    Parameters
    ----------
    px, py, pz : np.ndarray
        The positions of the SPH particles
    sigma : np.ndarray
        The sigma of the SPH particles (assuming a Gaussian Kernel)
    mass : np.ndarray
        The mass of the SPH particles
    mode : { 'fast', 'exact' }
        Whether to properly compute the integral over the kernel ('exact') or
        simply count the particles in each cell ('fast')

    Returns
    -------
    density : `~numpy.ndarray`
        The density in each cell.
    """

    # Find limits for each cell
    xmin, xmax, ymin, ymax, zmin, zmax = octree.limits

    if mode == 'exact':

        try:

            import multiprocessing as mp
            p = mp.Pool()
            p_map = p.map

            # Find number of processes that multiprocessing will use
            N = mp.cpu_count()

        except:

            p_map = map
            N = 1

        # Define indices
        idx = np.indices(xmin.shape)[0]

        # Split them up for multiprocessing
        size = int(np.ceil(len(idx) / float(N)))
        idx_split = [idx[i:i+size] for i in range(0, len(idx), size)]
        assert np.all(np.hstack(idx_split) == idx)

        # Construct tuple to send to multiprocessing
        arguments = []
        for idx_subset in idx_split:
            arguments.append((xmin[idx_subset], xmax[idx_subset],
                              ymin[idx_subset], ymax[idx_subset],
                              zmin[idx_subset], zmax[idx_subset],
                              px, py, pz,
                              sigma, mass))

        densities = p_map(discretize_wrapper, arguments)

        density = np.hstack(densities)

    elif mode == 'fast':

        raise NotImplementedError("")

        # need to take into account masses
        # optimize?
        # density = np.array([np.sum(x[4]) for x in particles])

    else:

        raise ValueError("Unknown mode: {0} - should be one of 'exact' or 'fast'".format(mode))

    # Reset density to zero in cells that are sub-divided. Here we have to use
    # where because otherwise if refined is an integer array it will not do the
    # right thing (it will only affect the first two density elements).
    density[np.where(octree.refined)] = 0.

    # Normalize by volume
    density = density / (xmax - xmin) / (ymax - ymin) / (zmax - zmin)

    return density
