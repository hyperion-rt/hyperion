# We now define a decorator for methods that needs access to the output HDF5
# file. This is necessary because h5py has issues with links pointing to
# groups that are in open files.

import h5py
from ..util.decorator import decorator


def on_the_fly_hdf5(f):
    return decorator(_on_the_fly_hdf5, f)


def _on_the_fly_hdf5(f, *args, **kwargs):
    preset = args[0].file is not None
    if not preset:
        args[0].file = h5py.File(args[0].filename, 'r')
    try:
        return f(*args, **kwargs)
    finally:
        if not preset:
            args[0].file.close()
            args[0].file = None
