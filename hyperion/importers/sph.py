import numpy as np


def refine(x, y, z, dx, dy, dz, px, py, pz, h):

    if len(px) < 3 and np.all(dx < h / 5.) and np.all(dy < h / 5.) and np.all(dz < h / 5.):
        return [False], [(px, py, pz)], [(x-dx, x+dx, y-dy, y+dy, z-dz, z+dz)]

    b_all = [True]
    p_all = [([],[],[])]
    l_all = [(x-dx, x+dx, y-dy, y+dy, z-dz, z+dz)]

    for xcomp in (np.less, np.greater):
        xsub = x - dx * 0.5 if xcomp is np.less else x + dx * 0.5

        for ycomp in (np.less, np.greater):
            ysub = y - dy * 0.5 if ycomp is np.less else y + dy * 0.5

            for zcomp in (np.less, np.greater):
                zsub = z - dz * 0.5 if zcomp is np.less else z + dz * 0.5

                keep = xcomp(px, x) & ycomp(py, y) & zcomp(pz, z)

                b, p, l = refine(xsub, ysub, zsub,
                                 dx * 0.5, dy * 0.5, dz * 0.5,
                                 px[keep], py[keep], pz[keep], h[keep])

                b_all += b
                p_all += p
                l_all += l

    return b_all, p_all, l_all


def construct_octree(x, y, z, dx, dy, dz, px, py, pz, h):
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
    h : np.ndarray
        The radii of the SPH particles

    Returns
    -------
    grid : `~hyperion.grid.octree_grid.OctreeGrid`
        The octree grid
    """

    from ..grid import OctreeGrid
    from ._discretize_sph import _discretize_sph_func

    refined, particles, limits = refine(x, y, z, dx, dy, dz, px, py, pz, h)

    octree = OctreeGrid(x, y, z, dx, dy, dz, refined)

    xmin, xmax, ymin, ymax, zmin, zmax = zip(*limits)

    sigmax, sigmay, sigmaz = h, h, h

    xmin = np.array(xmin)
    xmax = np.array(xmax)
    ymin = np.array(ymin)
    ymax = np.array(ymax)
    zmin = np.array(zmin)
    zmax = np.array(zmax)

    density = _discretize_sph_func(xmin, xmax, ymin, ymax, zmin, zmax, px, py, pz, sigmax, sigmay, sigmaz)

    octree['density'] = []
    octree['density'].append(density)

    return octree

# TODO: parallelize top-level search on 8 processes
