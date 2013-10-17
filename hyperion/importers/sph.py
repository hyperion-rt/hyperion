import numpy as np


def refine(x, y, z, dx, dy, dz, px, py, pz):

    if len(px) < 3:
        return [0], [(px, py, pz)], [(x-dx, x+dx, y-dy, y+dy, z-dz, z+dz)]

    b_all = [1]
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
                                 px[keep], py[keep], pz[keep])

                b_all += b
                p_all += p
                l_all += l

    return b_all, p_all, l_all


def construct_octree(x, y, z, dx, dy, dz, px, py, pz):
    """
    Construct an Octree grid from SPH particles

    Parameters
    ----------
    x, y, z : float
        The center of the top-level cell
    dx, dy, dz : float
        The half-width of the top-level cell
    px, py, pz : np.ndarray
        The SPH particles

    Returns
    -------
    grid : `~hyperion.grid.octree_grid.OctreeGrid`
        The octree grid
    """

    refined, particles, limits = refine(x, y, z, dx, dy, dz, px, py, pz)

    xmin, xmax, ymin, ymax, zmin, zmax = zip(*limits)
    
    print xmin


if __name__ == '__main__':

    np.random.seed(12345)

    N = 100

    XMID = 0.
    DX = 0.5
    YMID = 0.
    DY = 0.5
    ZMID = 0.
    DZ = 0.5
    
    px = np.random.uniform(XMID - DX, XMID + DX, N)
    py = np.random.uniform(YMID - DY, YMID + DY, N)
    pz = np.random.uniform(ZMID - DZ, ZMID + DZ, N)
    
    construct_octree(XMID, YMID, ZMID, DX, DY, DZ, px, py, pz)
    
