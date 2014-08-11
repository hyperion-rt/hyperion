from __future__ import print_function, division

import os

DATA = os.path.join(os.path.dirname(__file__), 'data')

# Helper functions.


# Check that all non-wall neighbours in voro are also in qhull. qhull
# could have extra neighbours due to the clipping of the domain
# by voro++.
def _neighbours_check(voro, qhull):
    return all([all([_ in set(filter(lambda n: n >= 0, tup[0]))] for _ in filter(lambda n: n >= 0, tup[1])) for tup in zip(qhull['neighbours'], voro['neighbours'])])

# Check that the volumes are positive and amount to the total domain volume.
def _volume_check(voro):
    c1 = all(v > 1E-13 for v in voro['volume'])
    c2 = abs(sum(voro['volume']) - 1.) < 1E-13
    return (c1 and c2)

# Check that the bounding boxes are consistent with the sites' positions.
def _bb_check(voro):
    retval = []
    for r in voro:
        c,m,M = r['coordinates'],r['bb_min'],r['bb_max']
        retval.append(all(c[i] >= m[i] and c[i] <= M[i] for i in range(0,3)))
    return all(retval)

# Test for internal consistency and consistency wrt
# the qhull version of the code. 
def test_consistency():
    from astropy.table import Table
    import random
    import numpy as np
    from .. import voronoi_helpers as vh
    # Three pickled neighbours tables generated with qhull.
    l = ['qhull_00.hdf5', 'qhull_01.hdf5', 'qhull_02.hdf5']
    for s in l:
        # The three qhull tables have been generated with seeds (0,1,2)
        # and with (10,100,1000) sites respectively.
        idx = l.index(s)
        random.seed(idx)
        sites_arr = np.array([(random.uniform(0, 1), random.uniform(
            0, 1), random.uniform(0, 1)) for i in range(0, 10 ** (idx + 1))])
        vg = vh.voronoi_grid(sites_arr, np.array([[0, 1.], [0, 1], [0, 1]]))
        voro = vg.neighbours_table
        qhull = Table.read(os.path.join(DATA, s),path='data')
        assert _neighbours_check(voro, qhull)
        assert _volume_check(voro)
        assert _bb_check(voro)
