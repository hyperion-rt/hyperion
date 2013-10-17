import numpy as np


def refine(xmin, xmax, ymin, ymax, zmin, zmax, px, py, pz):

    xmid = 0.5 * (xmin + xmax)
    ymid = 0.5 * (ymin + ymax)
    zmid = 0.5 * (zmin + zmax)

    if len(px) < 3:
        return [0], [(px, py, pz)], [(xmin, xmax, ymin, ymax, zmin, zmax)]

    b_all = [1]
    p_all = [([],[],[])]
    l_all = [(xmin, xmax, ymin, ymax, zmin, zmax)]

    k = (px < xmid) & (py < ymid) & (pz < zmid)
    b, p, l = refine(xmin, xmid, ymin, ymid, zmin, zmid, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    k = (px > xmid) & (py < ymid) & (pz < zmid)
    b, p, l = refine(xmid, xmax, ymin, ymid, zmin, zmid, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    k = (px < xmid) & (py > ymid) & (pz < zmid)
    b, p, l = refine(xmin, xmid, ymid, ymax, zmin, zmid, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    k = (px > xmid) & (py > ymid) & (pz < zmid)
    b, p, l = refine(xmid, xmax, ymid, ymax, zmin, zmid, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    k = (px < xmid) & (py < ymid) & (pz > zmid)
    b, p, l = refine(xmin, xmid, ymin, ymid, zmid, zmax, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    k = (px > xmid) & (py < ymid) & (pz > zmid)
    b, p, l = refine(xmid, xmax, ymin, ymid, zmid, zmax, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    k = (px < xmid) & (py > ymid) & (pz > zmid)
    b, p, l = refine(xmin, xmid, ymid, ymax, zmid, zmax, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    k = (px > xmid) & (py > ymid) & (pz > zmid)
    b, p, l = refine(xmid, xmax, ymid, ymax, zmid, zmax, px[k], py[k], pz[k])
    b_all += b
    p_all += p
    l_all += l

    return b_all, p_all, l_all


def construct_octree(xmin, xmax, ymin, ymax, zmin, zmax, px, py, pz):

    refine, particles, limits = refine(XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, x, y, z)

    xmin, xmax, ymin, ymax, zmin, zmax = zip(*limits)


if __name__ == '__main__':

    np.random.seed(12345)

    N = 100

    XMIN = 0.
    XMAX = 1.
    YMIN = 0.
    YMAX = 1.
    ZMIN = 0.
    ZMAX = 1.

    x = np.random.uniform(XMIN, XMAX, N)
    y = np.random.uniform(YMIN, YMAX, N)
    z = np.random.uniform(ZMIN, ZMAX, N)
