import h5py

from ..util.otf_hdf5 import on_the_fly_hdf5


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
