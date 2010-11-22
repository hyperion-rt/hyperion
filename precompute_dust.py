import glob
import os

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

np.seterr(all='ignore')

from hyperion.util.functions import filename2hdf5
from hyperion.dust import SimpleSphericalDust, prepare_emiss

for dustfile in glob.glob(os.path.join('dustfiles/', '*')):
    print "Processing %s" % dustfile
    dustfile_out = filename2hdf5(dustfile)
    d = SimpleSphericalDust(dustfile)
    d._extrapolate(1.e-3, 1.e5)
    d.write('hyperion/data/' + os.path.basename(dustfile_out))
    d.plot('hyperion/data/' + os.path.basename(dustfile_out).replace('.hdf5', '.png'))
