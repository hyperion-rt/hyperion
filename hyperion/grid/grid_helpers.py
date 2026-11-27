
import h5py
import numpy as np

from ..util.functions import is_numpy_array


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

        elif data.ndim == ndim + 2:
            # Frequency-resolved quantity, stored as (n_freq, n_pop, *grid)
            n_pop = data.shape[1]
            shape = data.shape[-ndim:]

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
        elif len(shape) == ndim + 2:
            n_pop = shape[1]
            shape = shape[-ndim:]
        else:
            raise Exception("Unexpected number of dimensions: %i" % len(shape))
    else:
        raise ValueError("Data should be a list or a Numpy array or an "
                         "external HDF5 link")

    return n_pop, shape


class GridViewArrayMixin(object):
    """
    Mixin for the grid view classes that gives access to the viewed quantity
    through the ``array`` attribute, including assigning to it, which in-place
    operators such as ``*=`` do.
    """

    def __getattr__(self, attribute):
        if attribute == 'array':
            return self.quantities[self.viewed_quantity]
        else:
            return super().__getattr__(attribute)

    def __setattr__(self, attribute, value):
        if attribute == 'array':
            current = self.quantities[self.viewed_quantity]
            # In-place operators such as *= assign back the same object
            if value is current:
                return
            if type(current) is list:
                raise TypeError("Can only set the array of a view of a single grid, e.g. grid['density'][0].array")
            if not is_numpy_array(value):
                raise TypeError("array should be a Numpy array")
            if value.shape != current.shape:
                raise ValueError("Quantity arrays do not have the right "
                                 "dimensions: %s instead of %s"
                                 % (value.shape, current.shape))
            # Modify in place so that the change propagates to the parent grid
            current[...] = value
        else:
            super().__setattr__(attribute, value)
