import h5py


from ..util.decorator import decorator

def on_the_fly_hdf5(f):
    return decorator(_on_the_fly_hdf5, f)


def _on_the_fly_hdf5(f, *args, **kwargs):
    preset = args[0].file is not None
    if not preset:
        args[0].file = h5py.File(args[0].filename, 'r')
    results = f(*args, **kwargs)
    if not preset:
        args[0].file.close()
        args[0].file = None
    return results


class GridOnDisk(object):

    def __init__(self, filename, path='/'):

        self.filename = filename
        self.path = path
        self.file = None

    @on_the_fly_hdf5
    def __contains__(self, item):
        return item in self.file['Quantities']

    def __getitem__(self, item):
        return GridQuantity(self, item)

    @property
    def link(self):
        return h5py.ExternalLink(self.filename, self.path)


class GridQuantity(object):

    def __init__(self, grid, quantity):
        self.filename = grid.filename
        self.path = grid.path
        self.quantity = quantity
        self.file = None

    @property
    @on_the_fly_hdf5
    def n_dust(self):
        return self.file['Quantities'][self.quantity].shape[0]
