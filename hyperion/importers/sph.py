import random
import numpy as np


def DEFAULT_STOPPING_CRITERION(x, y, z, dx, dy, dz, px, py, pz, sigma):
    return len(px) <= 2


def refine(x, y, z, dx, dy, dz, px, py, pz, sigma, levels_remaining, stopping_criterion):

    if stopping_criterion(x, y, z, dx, dy, dz, px, py, pz, sigma) or levels_remaining == 0:
        return [False], [(px, py, pz)], [(x-dx, x+dx, y-dy, y+dy, z-dz, z+dz)]

    b_all = [True]
    p_all = [([],[],[])]
    l_all = [(x-dx, x+dx, y-dy, y+dy, z-dz, z+dz)]

    px_pos = px > x
    py_pos = py > y
    pz_pos = pz > z

    for zcomp, zsub in ((~pz_pos, z - dz * 0.5),(pz_pos, z + dz  *0.5)):
        for ycomp, ysub in ((~py_pos, y - dy * 0.5),(py_pos, y + dy  *0.5)):
            for xcomp, xsub in ((~px_pos, x - dx * 0.5),(px_pos, x + dx  *0.5)):

                keep = xcomp & ycomp & zcomp

                b, p, l = refine(xsub, ysub, zsub,
                                 dx * 0.5, dy * 0.5, dz * 0.5,
                                 px[keep], py[keep], pz[keep], sigma[keep],
                                 levels_remaining - 1,
                                 stopping_criterion)

                b_all += b
                p_all += p
                l_all += l


    return b_all, p_all, l_all

# def refine_paralle_wrapper(args):
#     return refine(*args)
#
# def refine_parallel(x, y, z, dx, dy, dz, px, py, pz, h, level):
#
#     b_all = [True]
#     p_all = [([],[],[])]
#     l_all = [(x-dx, x+dx, y-dy, y+dy, z-dz, z+dz)]
#
#     px_pos = px > x
#     py_pos = py > y
#     pz_pos = pz > z
#
#     arguments = []
#
#     for xcomp, xsub in ((~px_pos, x - dx * 0.5),(px_pos, x + dx  *0.5)):
#         for ycomp, ysub in ((~py_pos, y - dy * 0.5),(py_pos, y + dy  *0.5)):
#             for zcomp, zsub in ((~pz_pos, z - dz * 0.5),(pz_pos, z + dz  *0.5)):
#
#                 keep = xcomp & ycomp & zcomp
#
#                 arguments.append((xsub, ysub, zsub,
#                                   dx * 0.5, dy * 0.5, dz * 0.5,
#                                   px[keep], py[keep], pz[keep], h[keep], level + 1))
#
#
#     import multiprocessing as mp
#     p = mp.Pool(processes=8)
#     results = zip(*p.map(refine_paralle_wrapper, arguments))
#
#     b_all = b_all + reduce(list.__add__, results[0])
#     p_all = p_all + reduce(list.__add__, results[1])
#     l_all = l_all + reduce(list.__add__, results[2])
#
#     return b_all, p_all, l_all

def discretize_wrapper(args):
    from ._discretize_sph import _discretize_sph_func
    return _discretize_sph_func(*args)


def construct_octree(x, y, z, dx, dy, dz, px, py, pz, sigma, mass, n_levels=None, stopping_criterion=DEFAULT_STOPPING_CRITERION):
    """
    Construct an Octree grid from SPH particles

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

    from ..grid import OctreeGrid

    if n_levels is None:
        n_levels = np.inf

    refined, particles, limits = refine(x, y, z, dx, dy, dz, px, py, pz, sigma, n_levels, stopping_criterion)

    octree = OctreeGrid(x, y, z, dx, dy, dz, refined)

    xmin, xmax, ymin, ymax, zmin, zmax = zip(*limits)

    xmin = np.array(xmin)
    xmax = np.array(xmax)
    ymin = np.array(ymin)
    ymax = np.array(ymax)
    zmin = np.array(zmin)
    zmax = np.array(zmax)

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
        arguments.append((xmin[idx_subset], xmax[idx_subset], ymin[idx_subset], ymax[idx_subset], zmin[idx_subset], zmax[idx_subset],
                          px, py, pz,
                          sigma, mass))

    densities = p_map(discretize_wrapper, arguments)

    density = np.hstack(densities)

    # Reset density to zero in cells that are sub-divided. Here we have to use
    # where because otherwise if refined is an integer array it will not do the
    # right thing (it will only affect the first two density elements).
    density[np.where(refined)] = 0.

    # Normalize by volume
    density = density / (xmax - xmin) / (ymax - ymin) / (zmax - zmin)

    octree['density'] = []
    octree['density'].append(density)

    return octree
