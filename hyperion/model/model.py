from __future__ import print_function, division

import os
import subprocess
import multiprocessing
from copy import deepcopy

import h5py
import numpy as np

from ..version import __version__
from ..util.functions import delete_file
from ..grid import CartesianGrid, SphericalPolarGrid, CylindricalPolarGrid, OctreeGrid, AMRGrid
from ..sources import PointSource, SphericalSource, ExternalSphericalSource, ExternalBoxSource, MapSource, PlaneParallelSource
from ..conf import RunConf, PeeledImageConf, BinnedImageConf, OutputConf
from ..util.constants import c
from ..util.functions import FreezableClass, link_or_copy, is_numpy_array, bool2str
from ..dust import SphericalDust
from ..util.logger import logger
from ..util.validator import validate_scalar

from .helpers import find_last_iteration
from .model_output import ModelOutput


class Configuration(FreezableClass):

    def __init__(self):
        self.output = OutputConf()
        self._freeze()


class Model(FreezableClass, RunConf):

    def __init__(self, name=None):
        '''
        Create a Model instance

        Parameters
        ----------
        name : str
            The name of the model. This can be a descriptive name to identify
            the model. It is used by Model.write() to generate an input
            filename if none is specified.
        '''

        self.conf = Configuration()

        self.name = name

        self.reset_dust()
        self.reset_sources()
        self.reset_images()

        self._minimum_temperature = None
        self._minimum_specific_energy = None

        self.set_monochromatic(False)

        self.grid = None

        self.filename = None

        self.init_run_conf()

        self._freeze()

    def reset_dust(self):
        self.dust = None

    def reset_sources(self):
        self.sources = []

    def reset_images(self):
        self.binned_output = None
        self.peeled_output = []

    def set_monochromatic(self, monochromatic, frequencies=None, wavelengths=None):
        '''
        Set whether to do the radiation transfer at specific
        frequencies/wavelengths.

        Parameters
        ----------
        monochromatic : bool
            Whether to carry out radiation transfer at specific frequencies
            or wavelengths
        frequencies : iterable of floats, optional
            The frequencies to compute the radiation transfer for, in Hz
        wavelengths : iterable of floats, optional
            The wavelengths to compute the radiation transfer for, in microns

        If `monochromatic` is True, then one of `frequencies` or
        `wavelengths` is required
        '''

        self._monochromatic = monochromatic

        if self._monochromatic:
            if wavelengths is not None and frequencies is not None:
                raise Exception("Cannot specify both frequencies and wavelengths")
            elif wavelengths is not None:
                frequencies = c / (np.array(wavelengths) * 1.e-4)
            elif frequencies is not None:
                frequencies = np.array(frequencies)
            else:
                raise Exception("Need to specify frequencies or wavelengths")
        else:
            if wavelengths is not None or frequencies is not None:
                raise Exception("Cannot specify monochromatic frequencies or wavelengths if monochromatic=False")

        self._frequencies = frequencies

        for images in self.peeled_output:
            images._monochromatic = True
            if type(images.wav_min) != int or type(images.wav_max) != int:
                images.set_wavelength_range(len(frequencies), 1, len(frequencies))
        if self.binned_output is not None:
            raise Exception("Binned images cannot be computed in monochromatic mode")

    def _write_monochromatic(self, group, compression=True, dtype=np.float64):
        group.attrs['monochromatic'] = bool2str(self._monochromatic)
        if self._monochromatic:
            group.create_dataset('frequencies', data=np.array(list(zip(self._frequencies)), dtype=[('nu', dtype)]), compression=compression)

    def use_geometry(self, filename):
        '''
        Use the grid from an existing output or input file

        Parameters
        ----------
        filename : str
            The file to read the grid from. This can be either the input or
            output file from a radiation transfer run.
        '''

        # Open existing file
        f = h5py.File(filename, 'r')

        # Get group pointing to grid
        if 'Grid' in f:
            g_grid = f['Grid']
        elif 'Grid' in f['Input']:
            if f['Input'].file != f.file:
                # Workaround for h5py bug - can't access link directly,
                # need to use file attribute
                g_grid = f['Input'].file[f['Input'].name]['Grid']
            else:
                g_grid = f['Input/Grid']
        else:
            raise Exception("No grid found in file: %s" % filename)

        # Determine grid type
        if g_grid['Geometry'].attrs['grid_type'].decode('utf-8') == 'car':
            grid = CartesianGrid()
        elif g_grid['Geometry'].attrs['grid_type'].decode('utf-8') == 'sph_pol':
            grid = SphericalPolarGrid()
        elif g_grid['Geometry'].attrs['grid_type'].decode('utf-8') == 'cyl_pol':
            grid = CylindricalPolarGrid()
        elif g_grid['Geometry'].attrs['grid_type'].decode('utf-8') == 'amr':
            grid = AMRGrid()
        elif g_grid['Geometry'].attrs['grid_type'].decode('utf-8') == 'oct':
            grid = OctreeGrid()
        else:
            raise NotImplemented("Cannot read geometry type %s" % g_grid['Geometry'].attrs['grid_type'].decode('utf-8'))

        # Read in the grid
        grid.read(g_grid, quantities=[])

        # Set the grid
        self.set_grid(grid)

        # Close the file
        f.close()

    def use_quantities(self, filename, quantities=['density', 'specific_energy'],
                       use_minimum_specific_energy=True, use_dust=True):
        '''
        Use physical quantities from an existing output file

        Parameters
        ----------
        filename : str
            The file to read the quantities from. This should be the output
            file of a radiation transfer run.
        quantities : list
            Which physical quantities to read in. Can include 'density' and
            'specific_energy'.
        copy : bool
            Whether to copy the quantities into the new input file, or whether
            to just link to them.
        use_minimum_specific_energy : bool
            Whether to also use the minimum specific energy values from the
            file specified
        use_dust : bool
            Whether to also use the dust properties from the file specified
        '''

        # Open existing file
        f = h5py.File(filename, 'r')

        # Find last iteration
        max_iteration = find_last_iteration(f)

        if max_iteration == 0:
            raise ValueError("No iterations found in file: %s" % filename)

        logger.info("Retrieving quantities from iteration %i of %s" % (max_iteration, filename))

        # Find path to file for link. For now, use absolute links.
        # In Model.write() we can always replace the links with
        # relative links if desired.
        file_path = os.path.abspath(filename)

        # Loop over quantities
        for quantity in ['density', 'specific_energy']:

            if quantity in quantities:

                # Set the path to the quantity
                if quantity in ['density']:
                    array_path = '/Input/Grid/Quantities/%s' % quantity
                else:
                    array_path = '/iteration_%05i/specific_energy' % max_iteration

                logger.info("Using %s from %s" % (quantity, filename))

                # Add quantity to grid
                self.grid[quantity] = h5py.ExternalLink(file_path, array_path)

        # Minimum specific energy
        if use_minimum_specific_energy:
            logger.info("Using minimum_specific_energy from %s" % filename)
            self.set_minimum_specific_energy([float(x) for x in f['/Input/Grid/Quantities'].attrs['minimum_specific_energy']])

        # Dust properties
        if use_dust:
            logger.info("Using dust properties from %s" % filename)
            self.dust = h5py.ExternalLink(file_path, '/Input/Dust')

        # Close the file
        f.close()

    def write(self, filename=None, compression=True, copy=True,
              absolute_paths=False, wall_dtype=float,
              physics_dtype=float, overwrite=True):
        '''
        Write the model input parameters to an HDF5 file

        Parameters
        ----------
        filename : str
            The name of the input file to write. If no name is specified, the
            filename is constructed from the model name.
        compression : bool
            Whether to compress the datasets inside the HDF5 file.
        copy : bool
            Whether to copy all external content into the input file, or
            whether to just link to external content.
        absolute_paths : bool
            If copy=False, then if absolute_paths is True, absolute filenames
            are used in the link, otherwise the path relative to the input
            file is used.
        wall_dtype : type
            Numerical type to use for wall positions.
        physics_dtype : type
            Numerical type to use for physical grids.
        overwrite : bool
            Whether to overwrite any pre-existing file
        '''

        # If no filename has been specified, use the model name to construct
        # one. If neither have been specified, raise an exception.
        if filename is None:
            if self.name is not None:
                filename = self.name + '.rtin'
            else:
                raise ValueError("filename= has not been specified and model "
                                 "has no name")

        # Remove previous file if it exists
        if overwrite and os.path.exists(filename):
            os.remove(filename)

        # Check that grid has been set up
        if self.grid is None:
            raise Exception("No coordinate grid has been set up")

        # Check that containing directory exists to give a more understandable
        # message than 'File does not exist', since this might confuse users
        # (it is the output directory that does not exist)
        if not os.path.dirname(filename) == "":
            if not os.path.exists(os.path.dirname(filename)):
                raise IOError("Directory %s does not exist" % \
                              os.path.dirname(filename))

        # Create output file
        delete_file(filename)
        root = h5py.File(filename, 'w')

        # Add Python version
        root.attrs['python_version'] = np.string_(__version__.encode('utf-8'))

        # Create all the necessary groups and sub-groups
        g_grid = root.create_group('Grid')
        g_sources = root.create_group('Sources')
        g_output = root.create_group('Output')
        g_peeled = g_output.create_group('Peeled')
        g_binned = g_output.create_group('Binned')

        # Output sources
        for i, source in enumerate(self.sources):
            if isinstance(source, MapSource):
                source.write(g_sources, 'source_%05i' % (i + 1), self.grid,
                             compression=compression, map_dtype=physics_dtype)
            else:
                source.write(g_sources, 'source_%05i' % (i + 1))

        # Output configuration for peeled images/SEDs
        for i, peel in enumerate(self.peeled_output):
            if not self._frequencies is None:
                if not peel._monochromatic:
                    raise Exception("Peeled images need to be set to monochromatic mode")
            peel.write(g_peeled.create_group('group_%05i' % (i + 1)))

        # Output configuration for binned images/SEDs
        if self.binned_output is not None:
            self.binned_output.write(g_binned.create_group('group_00001'))

        # Write monochromatic configuration
        self._write_monochromatic(root, compression=compression)

        # Write run-time and output configuration
        self.write_run_conf(root)
        self.conf.output.write(g_output)

        # Check self-consistency of grid
        self.grid._check_array_dimensions()

        # Write the geometry and physical quantity arrays to the input file
        self.grid.write(g_grid, copy=copy, absolute_paths=absolute_paths, compression=compression, physics_dtype=physics_dtype)

        if 'density' in self.grid:

            # Check if dust types are specified for each
            if self.dust is None:
                raise Exception("No dust properties specified")

            if isinstance(self.dust, h5py.ExternalLink):

                link_or_copy(root, 'Dust', self.dust, copy, absolute_paths=absolute_paths)

            elif isinstance(self.dust, h5py.Group):

                root.copy(self.dust, 'Dust')

            elif type(self.dust) == list:

                g_dust = root.create_group('Dust')

                if self.grid['density'].n_dust != len(self.dust):
                    raise Exception("Number of density grids should match number of dust types")

                # Output dust file, avoiding writing the same dust file multiple times
                present = {}
                for i, dust in enumerate(self.dust):

                    short_name = 'dust_%03i' % (i + 1)

                    if copy:

                        if isinstance(dust, basestring):
                            dust = SphericalDust(dust)

                        if dust.hash() in present:
                            g_dust[short_name] = h5py.SoftLink(present[dust.hash()])
                        else:
                            dust.write(g_dust.create_group(short_name))
                            present[dust.hash()] = short_name

                    else:

                        if type(dust) != str:
                            if dust.filename is None:
                                raise ValueError("Dust properties are not located in a file, so cannot link. Use copy=True or write the dust properties to a file first")
                            else:
                                dust = dust.filename

                        if absolute_paths:
                            path = os.path.abspath(dust)
                        else:
                            # Relative path should be relative to input file, not current directory.
                            path = os.path.relpath(dust, os.path.dirname(filename))

                        g_dust[short_name] = h5py.ExternalLink(path, '/')

            else:
                raise ValueError("Unknown type for dust attribute: %s" % type(self.dust))

            _n_dust = len(root['Dust'])

            # Write minimum specific energy
            if self._minimum_temperature is not None:

                if np.isscalar(self._minimum_temperature):
                    _minimum_temperature = [self._minimum_temperature for i in range(_n_dust)]
                elif len(self._minimum_temperature) != _n_dust:
                    raise Exception("Number of minimum_temperature values should match number of dust types")
                else:
                    _minimum_temperature = self._minimum_temperature

                _minimum_specific_energy = []
                for i, dust in enumerate(root['Dust']):
                    mo = SphericalDust(root['Dust'][dust]).mean_opacities
                    _minimum_specific_energy.append(mo._temperature2specific_energy(_minimum_temperature[i]))

            elif self._minimum_specific_energy is not None:

                if np.isscalar(self._minimum_specific_energy):
                    _minimum_specific_energy = [self._minimum_specific_energy for i in range(_n_dust)]
                elif len(self._minimum_specific_energy) != _n_dust:
                    raise Exception("Number of minimum_specific_energy values should match number of dust types")
                else:
                    _minimum_specific_energy = self._minimum_specific_energy

            else:

                _minimum_specific_energy = [0. for i in range(_n_dust)]

            g_grid['Quantities'].attrs["minimum_specific_energy"] = [float(x) for x in _minimum_specific_energy]

        else:

            root.create_group('Dust')

        root.close()

        self.filename = filename

    def add_point_source(self, *args, **kwargs):
        source = PointSource(*args, **kwargs)
        self.add_source(source)
        return source

    def add_spherical_source(self, *args, **kwargs):
        source = SphericalSource(*args, **kwargs)
        self.add_source(source)
        return source

    def add_external_spherical_source(self, *args, **kwargs):
        source = ExternalSphericalSource(*args, **kwargs)
        self.add_source(source)
        return source

    def add_external_box_source(self, *args, **kwargs):
        source = ExternalBoxSource(*args, **kwargs)
        self.add_source(source)
        return source

    def add_map_source(self, *args, **kwargs):
        source = MapSource(*args, **kwargs)
        self.add_source(source)
        return source

    def add_plane_parallel_source(self, *args, **kwargs):
        source = PlaneParallelSource(*args, **kwargs)
        self.add_source(source)
        return source

    def add_source(self, source):
        self.sources.append(source)

    def add_density_grid(self, density, dust, specific_energy=None, merge_if_possible=False):
        '''
        Add a density grid to the model

        Parameters
        ----------
        density : np.ndarray or grid quantity
            The density of the dust. This can be specified either as a 3-d
            Numpy array for cartesian, cylindrical polar, or spherical polar
            grids, as a 1-d array for octree grids, or as a grid quantity
            object for all grid types. Grid quantity objects are obtained by
            taking an instance of a grid class (e.g. ``AMRGrid``,
            ``CartesianGrid``, ...) and specifying the quantity as an index,
            e.g. ``amr['density']`` where ``amr`` is an ``AMRGrid`` object.
        dust : str or dust instance
            The dust properties, specified either as a string giving the
            filename of the dust properties, or a as an instance of a dust
            class (e.g. ``SphericalDust``, ``IsotropicDust``, ...).
        specific_energy : np.ndarray or grid quantity, optional
            The specific energy of the density grid. Note that in order for
            this to be useful, the number of initial iterations should be set
            to zero, otherwise these values will be overwritten after the
            first initial iteration.
        merge_if_possible : bool
            Whether to merge density arrays that have the same dust type
        '''

        # Check that grid has been previously defined
        if not self.grid:
            raise Exception("A coordinate system/grid has to be defined before adding a density grid")

        # Check whether grid geometries are the same
        self.grid._check_array_dimensions(density)
        if specific_energy is not None:
            self.grid._check_array_dimensions(specific_energy)

        # Check whether all densities are zero
        if np.all(density == 0.):
            logger.info("All density values are zero - ignoring density grid")
            return

        # Check consistency between density list size and specific energy list size
        if 'density' in self.grid:
            if specific_energy is not None and 'specific_energy' not in self.grid:
                raise Exception("Cannot add specific energy as it was not added for previous density arrays")
        else:
            self.dust = []
            self.grid['density'] = []
            if specific_energy is not None:
                self.grid['specific_energy'] = []

        # Check whether the density can be added to an existing one
        if merge_if_possible:

            # Only consider this if the specific energy is not specified
            if specific_energy is None:

                if isinstance(dust, basestring):

                    if dust in self.dust:

                        logger.info("Merging densities (identical filenames)")

                        ip = self.dust.index(dust)
                        self.grid['density'][ip].add(density)
                        return

                else:

                    dust_hashes = []
                    for d in self.dust:
                        if not isinstance(d, basestring):
                            dust_hashes.append(d.hash())
                        else:
                            dust_hashes.append(None)

                    if dust.hash() in dust_hashes:

                        logger.info("Merging densities (identical hashes)")

                        ip = dust_hashes.index(dust.hash())
                        self.grid['density'][ip].add(density)
                        return

        # Set the density and dust
        self.grid['density'].append(density)
        self.dust.append(dust)

        # Set specific energy if specified
        if specific_energy is not None:
            self.grid['specific_energy'].append(specific_energy)

    def set_cartesian_grid(self, x_wall, y_wall, z_wall):
        self.set_grid(CartesianGrid(x_wall, y_wall, z_wall))

    def set_cylindrical_polar_grid(self, w_wall, z_wall, p_wall):
        self.set_grid(CylindricalPolarGrid(w_wall, z_wall, p_wall))

    def set_spherical_polar_grid(self, r_wall, t_wall, p_wall):
        self.set_grid(SphericalPolarGrid(r_wall, t_wall, p_wall))

    def set_octree_grid(self, x, y, z, dx, dy, dz, refined):
        self.set_grid(OctreeGrid(x, y, z, dx, dy, dz, refined))

    def set_amr_grid(self, description):
        self.set_grid(AMRGrid(description))

    def set_grid(self, grid):
        if isinstance(grid, AMRGrid):
            self.grid = AMRGrid(grid)
        else:
            self.grid = deepcopy(grid)

    def add_peeled_images(self, **kwargs):
        self.peeled_output.append(PeeledImageConf(**kwargs))
        return self.peeled_output[-1]

    def add_binned_images(self, **kwargs):
        if self.binned_output:
            raise Exception("Only one set of binned images can be set at this time")
        else:
            self.binned_output = BinnedImageConf(**kwargs)
            return self.binned_output

    def set_minimum_temperature(self, temperature):
        '''
        Set the minimum temperature for the dust

        Parameters
        ----------
        temperature : float, list, tuple, or Numpy array

        Notes
        -----
        This method should not be used in conjunction with
        ``set_minimum_specific_energy`` - only one of the two should be used.
        '''
        if np.isscalar(temperature):
            validate_scalar('temperature', temperature, domain='positive')
        elif type(temperature) in [list, tuple] or is_numpy_array(temperature):
            temperature = list(temperature)
            for value in temperature:
                validate_scalar('temperature', value, domain='positive')
        if self._minimum_specific_energy is not None:
            raise Exception("minimum specific_energy has already been set")
        self._minimum_temperature = temperature

    def set_minimum_specific_energy(self, specific_energy):
        '''
        Set the minimum specific energy for the dust

        Parameters
        ----------
        specific_energy : float, list, tuple, or Numpy array

        Notes
        -----
        This method should not be used in conjunction with
        ``set_minimum_temperature`` - only one of the two should be used.
        '''
        if np.isscalar(specific_energy):
            validate_scalar('specific_energy', specific_energy, domain='positive')
        elif type(specific_energy) in [list, tuple] or is_numpy_array(specific_energy):
            specific_energy = list(specific_energy)
            for value in specific_energy:
                validate_scalar('specific_energy', value, domain='positive')
        if self._minimum_temperature is not None:
            raise Exception("minimum temperature has already been set")
        self._minimum_specific_energy = specific_energy

    def run(self, filename=None, logfile=None, mpi=False, n_processes=multiprocessing.cpu_count(), overwrite=False):

        if self.filename is None:
            raise ValueError("Input file does not exist - write() needs to be called before run()")

        if mpi:
            option = '-m {0}'.format(n_processes)
        else:
            option = ''

        input_file = self.filename
        if filename is None:
            if '.rtin' in self.filename:
                output_file = self.filename.replace('.rtin', '.rtout')
            else:
                output_file = self.filename + '.rtout'
        else:
            output_file = filename

        if overwrite and os.path.exists(output_file):
            os.remove(output_file)

        if logfile:
            flog = open(logfile, 'wb')
            returncode = subprocess.call('hyperion %s %s %s' % (option, input_file, output_file), stdout=flog, stderr=flog, shell=True)
        else:
            returncode = subprocess.call('hyperion %s %s %s' % (option, input_file, output_file), shell=True)

        if returncode != 0:
            raise SystemExit("An error occurred, and the run did not complete")

        # Return handle to output file
        return ModelOutput(output_file)
