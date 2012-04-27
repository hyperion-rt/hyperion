from __future__ import print_function, division

import struct
import hashlib

import h5py
import numpy as np

from ..util.functions import FreezableClass, is_numpy_array, link_or_copy
from ..util.logger import logger
from .grid_helpers import single_grid_dims


class OctreeGrid(FreezableClass):

    def __init__(self, *args):

        self.shape = None

        self.x = None
        self.y = None
        self.z = None

        self.dx = None
        self.dy = None
        self.dz = None

        self.refined = None

        self.quantities = {}

        self._freeze()

        if len(args) > 0:
            self.set_walls(*args)

    def set_walls(self, x, y, z, dx, dy, dz, refined):

        if not np.isscalar(x):
            raise ValueError("x should be a scalar value")
        if not np.isreal(x):
            raise ValueError("x should be a numerical value")

        if not np.isscalar(y):
            raise ValueError("y should be a scalar value")
        if not np.isreal(y):
            raise ValueError("y should be a numerical value")

        if not np.isscalar(z):
            raise ValueError("z should be a scalar value")
        if not np.isreal(z):
            raise ValueError("z should be a numerical value")

        if not np.isscalar(dx):
            raise ValueError("dx should be a scalar value")
        if not np.isreal(dx):
            raise ValueError("dx should be a numerical value")

        if not np.isscalar(dy):
            raise ValueError("dy should be a scalar value")
        if not np.isreal(dy):
            raise ValueError("dy should be a numerical value")

        if not np.isscalar(dz):
            raise ValueError("dz should be a scalar value")
        if not np.isreal(dz):
            raise ValueError("dz should be a numerical value")

        if type(refined) in [list, tuple]:
            refined = np.array(refined)

        if not is_numpy_array(refined) or refined.ndim != 1 or refined.dtype != bool:
            raise ValueError("refined should be a 1-D boolean sequence")

        self.x = x
        self.y = y
        self.z = z

        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.refined = refined

        self.shape = (1, 1, len(refined))

    def __getattr__(self, attribute):
        if attribute == 'n_dust':
            n_dust = None
            for quantity in self.quantities:
                n_dust_q, shape_q = single_grid_dims(self.quantities[quantity])
                if n_dust is None:
                    n_dust = n_dust_q
                else:
                    if n_dust != n_dust_q:
                        raise ValueError("Not all dust lists in the grid have the same size")
            return n_dust
        else:
            return FreezableClass.__getattribute__(self, attribute)

    def _check_array_dimensions(self, array=None):
        '''
        Check that a grid's array dimensions agree with this grid's metadata

        Parameters
        ----------
        array: np.ndarray or list of np.ndarray, optional
            The array for which to test the dimensions. If this is not
            specified, this method performs a self-consistency check of array
            dimensions and meta-data.
        '''

        n_pop_ref = None

        for quantity in self.quantities:

            n_pop, shape = single_grid_dims(self.quantities[quantity])

            if shape != self.shape:
                raise ValueError("Quantity arrays do not have the right "
                                 "dimensions: %s instead of %s"
                                 % (shape, self.shape))

            if n_pop is not None:
                if n_pop_ref is None:
                    n_pop_ref = n_pop
                elif n_pop != n_pop_ref:
                    raise ValueError("Not all dust lists in the grid have the same size")

    def read(self, group, quantities='all'):
        '''
        Read in an oct-tree grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to read the grid from
        quantities: 'all' or list
            Which physical quantities to read in. Use 'all' to read in all
            quantities or a list of strings to read only specific quantities.
        '''

        # Extract HDF5 groups for geometry and physics

        g_geometry = group['Geometry']
        g_quantities = group['Quantities']

        # Read in geometry

        if g_geometry.attrs['grid_type'].decode('utf-8') != 'oct':
            raise ValueError("Grid is not an oct-tree")

        self.set_walls(g_geometry.attrs['x'],
                       g_geometry.attrs['y'],
                       g_geometry.attrs['z'],
                       g_geometry.attrs['dx'],
                       g_geometry.attrs['dy'],
                       g_geometry.attrs['dz'],
                       g_geometry['cells']['refined'].astype(bool))

        # Read in physical quantities
        if quantities is not None:
            for quantity in g_quantities:
                if quantities == 'all' or quantity in quantities:
                    array = np.array(g_quantities[quantity])
                    if array.ndim == 4:  # if array is 4D, it is a list of 3D arrays
                        self.quantities[quantity] = [array[i] for i in range(array.shape[0])]
                    else:
                        self.quantities[quantity] = array

        # Check that advertised hash matches real hash
        if g_geometry.attrs['geometry'].decode('utf-8') != self.get_geometry_id():
            raise Exception("Calculated geometry hash does not match hash in file")

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

    def write(self, group, quantities='all', copy=True, absolute_paths=False, compression=True, wall_dtype=float, physics_dtype=float):
        '''
        Write out the oct-tree grid

        Parameters
        ----------
        group: h5py.Group
            The HDF5 group to write the grid to
        quantities: 'all' or list
            Which physical quantities to write out. Use 'all' to write out all
            quantities or a list of strings to write only specific quantities.
        compression: bool
            Whether to compress the arrays in the HDF5 file
        wall_dtype: type
            The datatype to use to write the wall positions (ignored for this kind of grid)
        physics_dtype: type
            The datatype to use to write the physical quantities
        '''

        # Create HDF5 groups if needed

        if 'Geometry' not in group:
            g_geometry = group.create_group('Geometry')
        else:
            g_geometry = group['Geometry']

        if 'Quantities' not in group:
            g_quantities = group.create_group('Quantities')
        else:
            g_quantities = group['Quantities']

        # Write out geometry

        g_geometry.attrs['grid_type'] = np.string_('oct'.encode('utf-8'))
        g_geometry.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

        g_geometry.attrs['x'] = self.x
        g_geometry.attrs['y'] = self.y
        g_geometry.attrs['z'] = self.z
        g_geometry.attrs['dx'] = self.dx
        g_geometry.attrs['dy'] = self.dy
        g_geometry.attrs['dz'] = self.dz

        dset = g_geometry.create_dataset("cells", data=np.array(list(zip(self.refined)), dtype=[('refined', np.int32)]), compression=compression)

        # Self-consistently check geometry and physical quantities
        self._check_array_dimensions()

        # Write out physical quantities

        for quantity in self.quantities:
            if quantities == 'all' or quantity in quantities:
                if isinstance(self.quantities[quantity], h5py.ExternalLink):
                    link_or_copy(g_quantities, quantity, self.quantities[quantity], copy, absolute_paths=absolute_paths)
                else:
                    dset = g_quantities.create_dataset(quantity, data=self.quantities[quantity],
                                                       compression=compression,
                                                       dtype=physics_dtype)
                    dset.attrs['geometry'] = np.string_(self.get_geometry_id().encode('utf-8'))

    def get_geometry_id(self):
        geo_hash = hashlib.md5()
        geo_hash.update(struct.pack('>d', self.x))
        geo_hash.update(struct.pack('>d', self.y))
        geo_hash.update(struct.pack('>d', self.z))
        geo_hash.update(struct.pack('>d', self.dx))
        geo_hash.update(struct.pack('>d', self.dy))
        geo_hash.update(struct.pack('>d', self.dz))
        geo_hash.update(self.refined)
        return geo_hash.hexdigest()

    def __getitem__(self, item):
        return OctreeGridView(self, item)

    def __setitem__(self, item, value):
        if isinstance(value, OctreeGridView):
            if self.refined is None:
                logger.warn("No geometry in target grid - copying from original grid")
                self.set_walls(value.x, value.y, value.z, value.dx, value.dy, value.dz, value.refined)
            self.quantities[item] = value.quantities[value.viewed_quantity]
        elif isinstance(value, h5py.ExternalLink):
            self.quantities[item] = value
        elif value == []:
            self.quantities[item] = []
        else:
            raise ValueError('value should be an empty list, and ExternalLink, or a OctreeGridView instance')

    def __contains__(self, item):
        return self.quantities.__contains__(item)

    def reset_quantities(self):
        self.quantities = {}


class OctreeGridView(OctreeGrid):

    def __init__(self, grid, quantity):
        self.viewed_quantity = quantity
        OctreeGrid.__init__(self)
        self.set_walls(grid.x, grid.y, grid.z, grid.dx, grid.dy, grid.dz, grid.refined)
        self.quantities = {quantity: grid.quantities[quantity]}

    def append(self, grid):
        '''
        Used to append quantities from another grid

        Parameters
        ----------
        grid: 1D Numpy array or OctreeGridView instance
            The grid to copy the quantity from
        '''
        if isinstance(grid, OctreeGridView):
            self.quantities[self.viewed_quantity].append(grid.quantities[grid.viewed_quantity])
        elif type(grid) is np.ndarray:
            if grid.ndim == 1:
                grid = grid.reshape((1, 1, grid.shape[0]))
            self.quantities[self.viewed_quantity].append(grid)
        else:
            raise ValueError("grid should be a Numpy array or a OctreeGridView object")
