from __future__ import print_function, division

import h5py
import numpy as np


def single_grid_dims(data, ndim=3):
    '''
    Find the number of populations

    Parameters
    ----------
    data : list or tuple or np.ndarray or h5py.ExternalLink
        The data to find the number of populations and shape for

    Returns
    -------
    n_pop : int
        Number of (dust) populations
    shape : tuple
        The dimensions of the grid
    '''

    if type(data) in [list, tuple]:

        n_pop = len(data)
        shape = None
        for item in data:
            if shape is None:
                shape = item.shape
            elif item.shape != shape:
                raise ValueError("Grids in list/tuple should have the same "
                                 "dimensions")
        if shape is not None and len(shape) != ndim:
            raise ValueError("Grids should be %i-dimensional" % ndim)

    elif isinstance(data, np.ndarray):

        if data.ndim == ndim:
            n_pop = None
            shape = data.shape
        elif data.ndim == ndim + 1:
            n_pop = data.shape[0]
            shape = data[0].shape
        else:
            raise Exception("Unexpected number of dimensions: %i" % data.ndim)

    elif isinstance(data, h5py.ExternalLink):

        f = h5py.File(data.filename, 'r')
        shape = f[data.path].shape
        f.close()

        if len(shape) == ndim:
            n_pop = None
        elif len(shape) == ndim + 1:
            n_pop = shape[0]
            shape = shape[1:]
        else:
            raise Exception("Unexpected number of dimensions: %i" % len(shape))
    else:
        raise ValueError("Data should be a list or a Numpy array or an "
                         "external HDF5 link")

    return n_pop, shape
