from __future__ import print_function, division

import os
import glob
import hashlib

import h5py

from .version import __version__

data_dir = __path__[0] + '/data/'

datafiles = {}
for datafile in glob.glob(os.path.join(data_dir, '*.hdf5')):
    f = h5py.File(datafile)
    hash = f.attrs['asciimd5'].decode('utf-8')
    datafiles[hash] = os.path.abspath(datafile)

def get_HDF5_datafile(filename):

    h = hashlib.md5(file(filename,'rb').read()).hexdigest()

    if h in datafiles:
        return datafiles[h]
    else:
        raise Exception("File does not exist")
