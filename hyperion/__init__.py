import hyperion.densities
import hyperion.dust
import hyperion.util
import hyperion.atmos
import hyperion.model

__version__ = '0.6.4'

data_dir = __path__[0] + '/data/'

import os
import glob
import hashlib
import h5py

datafiles = {}
for datafile in glob.glob(os.path.join(data_dir, '*.hdf5')):
    f = h5py.File(datafile)
    hash = f.attrs['asciimd5']
    datafiles[hash] = os.path.abspath(datafile)
    
def get_HDF5_datafile(filename):
    
    h = hashlib.md5(file(filename,'rb').read()).hexdigest()

    if h in datafiles:
        return datafiles[h]
    else:
        raise Exception("File does not exist")
