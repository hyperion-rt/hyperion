import warnings

import h5py
import numpy as np


class NaNWarning(UserWarning):
    pass


def number_nan(array):
    if (np.isscalar(array) and np.isreal(array)) or array.dtype.kind in ['i', 'f']:
        return np.sum(np.isnan(array))
    else:
        return 0


def check_for_nans(handle):
    """
    Check for NaNs anywhere in an HDF5 group
    """

    # List datasets and groups
    datasets = []
    groups = []

    def func(name, obj):
        if isinstance(obj, h5py.Dataset):
            datasets.append(name)
        elif isinstance(obj, h5py.Group):
            groups.append(name)
    handle.visititems(func)

     # Visit order is not guaranteed, so sort
    datasets.sort()
    groups.sort()

    # Loop over datasets to check for NaN values
    for d in datasets:
        array = np.array(handle[d])
        if array.dtype.kind == 'V':  # structured array
            for col in array.dtype.fields:
                n_nan = number_nan(array[col])
                if n_nan > 0:
                    warnings.warn("{0} NaN value(s) encountered in field {1} of dataset '{2}'".format(n_nan, col, d), NaNWarning)
        else:
            n_nan = number_nan(array)
            if n_nan > 0:
                warnings.warn("{0} NaN value(s) encountered in dataset '{1}'".format(n_nan, d), NaNWarning)

    # Loop over all groups and datasets to check attributes
    for item in ['/'] + datasets + groups:

        # Find all attributes
        attributes = list(handle[item].attrs.keys())
        attributes.sort()

        for a in attributes:
            n_nan = number_nan(handle[item].attrs[a])
            if n_nan > 0:
                warnings.warn("{0} NaN value(s) encountered in attribute '{1}' of object '{2}'".format(n_nan, a, item), NaNWarning)
