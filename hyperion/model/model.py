from __future__ import print_function, division

import os
import subprocess
import multiprocessing
from copy import deepcopy

import h5py
import numpy as np

from ..version import __version__
from ..util.functions import delete_file
from ..grid import CartesianGrid, SphericalPolarGrid, CylindricalPolarGrid, OctreeGrid, AMRGrid, VoronoiGrid, GridOnDisk
from ..sources import PointSource, PointSourceCollection, SphericalSource, ExternalSphericalSource, ExternalBoxSource, MapSource, PlaneParallelSource, read_source
from ..conf import RunConf, PeeledImageConf, BinnedImageConf, OutputConf
from ..util.constants import c
from ..util.functions import FreezableClass, link_or_copy, is_numpy_array, bool2str, str2bool
from ..dust import SphericalDust
from astropy import log as logger
from ..util.validator import validate_scalar

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

        super(Model, self).__init__()

        self._freeze()

    def reset_dust(self):
        self.dust = None

    def reset_sources(self):
        self.sources = []

    def reset_images(self):
        self.binned_output = None
        self.peeled_output = []

    def set_monochromatic(self, monochromatic, wavelengths=None, frequencies=None):
        '''
        Set whether to do the radiation transfer at specific wavelengths.

        Parameters
        ----------
        monochromatic : bool
            Whether to carry out radiation transfer at specific frequencies
            or wavelengths
        wavelengths : iterable of floats, optional
            The wavelengths to compute the radiation transfer for, in microns

        If `monochromatic` is True then `wavelengths` is required
        '''

        if frequencies is not None:
            logger.warn("The frequencies= option will soon be deprecated - please specify wavelengths instead.")

        self._monochromatic = monochromatic

        if self._monochromatic:

            if wavelengths is not None:
                self._frequencies = c / (np.array(wavelengths) * 1.e-4)
            elif frequencies is not None:
                self._frequencies = frequencies
            else:
                raise Exception("Need to specify wavelengths")

            for images in self.peeled_output:
                images._set_monochromatic(True, frequencies=self._frequencies)

            if self.binned_output is not None:
                raise Exception("Binned images cannot be computed in monochromatic mode")

        else:

            if wavelengths is not None:
                raise Exception("Cannot specify monochromatic wavelengths if monochromatic=False")

            self._frequencies = None

    def _read_monochromatic(self, group):
        self._monochromatic = str2bool(group.attrs['monochromatic'])
        if self._monochromatic:
            self._frequencies = np.array(group['frequencies']['nu'])

    def _write_monochromatic(self, group, compression=True, dtype=np.float64):
        group.attrs['monochromatic'] = bool2str(self._monochromatic)
        if self._monochromatic:
            group.create_dataset('frequencies', data=np.array(list(zip(self._frequencies)), dtype=[('nu', dtype)]), compression=compression)

    @classmethod
    def read(cls, filename, only_initial=True):
        """
        Read in a previous model file

        This can be used to read in a previous input file, or the input in an
        output file (which is possible because the input to a model is stored
        or linked in an output file).

        If you are interested in re-using the final specific energy (and
        final density, if present) of a previously run model, you can use::

            >>> m = Model.read('previous_model.rtout', only_initial=False)

        Parameters
        ----------
        filename : str
            The name of the file to read the input data from
        only_initial : bool, optional
            Whether to use only the initial quantities, or whether the final
            specific energy (and optionally density) from the previous model
            can be used. By default, only the input density (and specific
            energy, if present) are read in.
        """

        self = cls()
        self.use_geometry(filename)
        self.use_quantities(filename, only_initial=only_initial)
        self.use_sources(filename)
        self.use_monochromatic_config(filename)
        self.use_run_config(filename)
        self.use_image_config(filename)
        self.use_output_config(filename)
        return self

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
        elif g_grid['Geometry'].attrs['grid_type'].decode('utf-8') == 'vor':
            grid = VoronoiGrid()
        else:
            raise NotImplemented("Cannot read geometry type %s" % g_grid['Geometry'].attrs['grid_type'].decode('utf-8'))

        # Read in the grid
        grid.read(g_grid, quantities=[])

        # Set the grid
        self.set_grid(grid)

        # Close the file
        f.close()

    def use_quantities(self, filename, quantities=['density', 'specific_energy'],
                       use_minimum_specific_energy=True, use_dust=True, copy=True,
                       only_initial=False):
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
        copy : bool
            Whether to read in a copy of the data. If set to False, then the
            physical quantities will only be links to the specified HDF5 file,
            and can therefore not be modified.
        only_initial : bool, optional
            Whether to use only the initial quantities, or whether the final
            specific energy (and optionally density) from the previous model
            can be used. By default, only the input density (and specific
            energy, if present) are read in.
        '''

        from .helpers import find_last_iteration

        # Open existing file
        f = h5py.File(filename, 'r')

        # Find path to file for link. For now, use absolute links.
        # In Model.write() we can always replace the links with
        # relative links if desired.
        file_path = os.path.abspath(filename)

        logger.info("Retrieving quantities from %s" % filename)

        quantities_path = {}

        if 'Input' in f:

            # Find last iteration
            max_iteration = find_last_iteration(f)

            if only_initial:
                logger.info("Reading input quantities")
            elif max_iteration == 0:
                logger.warn("No iterations found in file - only the input quantities will be used")
                last_iteration = None
            else:
                logger.info("Retrieving quantities from iteration %i" % max_iteration)
                last_iteration = 'iteration_{0:05d}'.format(max_iteration)

            if 'density' in quantities:
                if only_initial or last_iteration is None or 'density' not in f[last_iteration]:
                    quantities_path['density'] = '/Input/Grid/Quantities'
                else:
                    quantities_path['density'] = last_iteration

            if 'specific_energy' in quantities:
                if only_initial or last_iteration is None:
                    if 'specific_energy' in f['/Input/Grid/Quantities']:
                        quantities_path['specific_energy'] = '/Input/Grid/Quantities'
                else:
                    quantities_path['specific_energy'] = last_iteration

            # Minimum specific energy
            if use_minimum_specific_energy:
                minimum_specific_energy_path = '/Input/Grid/Quantities'

            if use_dust:
                dust_path = '/Input/Dust'

        else:

            if 'density' in quantities:
                quantities_path['density'] = '/Grid/Quantities'

            if 'specific_energy' in quantities:
                if 'specific_energy' in f['/Grid/Quantities']:
                    quantities_path['specific_energy'] = '/Grid/Quantities'

            # Minimum specific energy
            if use_minimum_specific_energy:
                minimum_specific_energy_path = '/Grid/Quantities'

            if use_dust:
                dust_path = '/Dust'

        # Now extract the quantities
        for quantity in quantities_path:

            logger.info("Using {quantity} from {path} in {filename}".format(quantity=quantity, path=quantities_path[quantity], filename=filename))

            # Add quantity to grid
            if copy:
                self.grid.read_quantities(f[quantities_path[quantity]], quantities=[quantity])
            else:
                self.grid[quantity] = h5py.ExternalLink(file_path, quantities_path[quantity] + '/' + quantity)

        # Minimum specific energy
        if use_minimum_specific_energy and 'minimum_specific_energy' in f[minimum_specific_energy_path].attrs:
            logger.info("Using minimum_specific_energy from {filename}".format(filename=filename))
            self.set_minimum_specific_energy([float(x) for x in f[minimum_specific_energy_path].attrs['minimum_specific_energy']])

        # Dust properties
        if use_dust:
            logger.info("Using dust properties from {filename}".format(filename=filename))
            if copy:
                self.dust = [SphericalDust(f[dust_path][name]) for name in f[dust_path]]
            else:
                self.dust = h5py.ExternalLink(file_path, dust_path)

        # Close the file
        f.close()

    def use_sources(self, filename):
        '''
        Use sources from an existing output file

        Parameters
        ----------
        filename : str
            The file to read the sources from. This should be the input or
            output file of a radiation transfer run.
        '''

        logger.info("Retrieving sources from %s" % filename)

        # Open existing file
        f = h5py.File(filename, 'r')

        # Get a pointer to the group with the sources
        if 'Sources' in f:
            g_sources = f['/Sources/']
        else:
            g_sources = f['/Input/Sources/']

        # Loop over sources
        for source in g_sources:
            self.add_source(read_source(g_sources[source]))

        # Close the file
        f.close()

    def use_monochromatic_config(self, filename):
        '''
        Use monochromatic configuration from an existing output file.

        Parameters
        ----------
        filename : str
            The file to read the configuration from. This should be the input or
            output file of a radiation transfer run.
        '''

        logger.info("Retrieving monochromatic configuration from %s" % filename)

        # Open existing file
        f = h5py.File(filename, 'r')

        # Get a pointer to the group with the sources
        if 'Input' in f:
            g = f['/Input']
        else:
            g = f

        # Read in monochromatic configuration
        self._read_monochromatic(g)

        # Close the file
        f.close()

    def use_run_config(self, filename):
        '''
        Use runtime configuration from an existing output or input file

        Parameters
        ----------
        filename : str
            The file to read the parameters from. This can be either the input
            or output file from a radiation transfer run.
        '''

        # need to do this here because n_photons will depend on monochromatic vs not
        self.use_monochromatic_config(filename)

        logger.info("Retrieving runtime configuration from %s" % filename)

        # Open existing file
        f = h5py.File(filename, 'r')

        # Get a pointer to the group with the sources
        if 'Input' in f:
            g_par = f['/Input/']
        else:
            g_par = f

        # Read in runtime configuration
        self.read_run_conf(g_par)

    def use_image_config(self, filename):
        '''
        Use image configuration from an existing output or input file

        Parameters
        ----------
        filename : str
            The file to read the parameters from. This can be either the input
            or output file from a radiation transfer run.
        '''

        # need to do this here because image wavelength interval will depend on monochromatic vs not
        self.use_monochromatic_config(filename)

        logger.info("Retrieving image configuration from %s" % filename)

        # Open existing file
        f = h5py.File(filename, 'r')

        # Get a pointer to the group with the sources
        if 'Output' in f:
            g_image = f['/Output/']
        else:
            g_image = f['/Input/Output/']

        # Read in binned images
        if 'n_theta' in g_image['Binned']:
            self.binned_output = BinnedImageConf.read(g_image['Binned'])

        # Read in peeled images
        for peeled in g_image['Peeled']:
            self.peeled_output.append(PeeledImageConf.read(g_image['Peeled'][peeled]))

    def use_output_config(self, filename):
        '''
        Use output configuration from an existing output or input file

        Parameters
        ----------
        filename : str
            The file to read the parameters from. This can be either the input
            or output file from a radiation transfer run.
        '''

        logger.info("Retrieving output configuration from %s" % filename)

        # Open existing file
        f = h5py.File(filename, 'r')

        # Get a pointer to the group with the sources
        if 'Output' in f:
            g_output = f['/Output/']
        else:
            g_output = f['/Input/Output/']

        # Read in output configuration
        self.conf.output.read(g_output)

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
                raise IOError("Directory %s does not exist" %
                              os.path.dirname(filename))

        # Create output file
        delete_file(filename)
        root = h5py.File(filename, 'w')

        # Add Python version
        root.attrs['python_version'] = np.string_(__version__.encode('utf-8'))

        # Create all the necessary groups and sub-groups
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

        if isinstance(self.grid, GridOnDisk):

            g_grid = link_or_copy(root, 'Grid', self.grid.link, copy=copy, absolute_paths=absolute_paths)

        else:

            # Create group
            g_grid = root.create_group('Grid')

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
                            if dust._file is None:
                                raise ValueError("Dust properties are not located in a file, so cannot link. Use copy=True or write the dust properties to a file first")
                            else:
                                # Check that has still matches file
                                if dust.hash() != dust._file[1]:
                                    raise ValueError("Dust properties have been modified since "
                                                     "being read in, so cannot link to dust file "
                                                     "on disk. You can solve this by writing out "
                                                     "the dust properties to a new file, or by "
                                                     "using copy=True.")
                                dust = dust._file[0]

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
                    d = SphericalDust(root['Dust'][dust])
                    _minimum_specific_energy.append(d.temperature2specific_energy(_minimum_temperature[i]))

            elif self._minimum_specific_energy is not None:

                if np.isscalar(self._minimum_specific_energy):
                    _minimum_specific_energy = [self._minimum_specific_energy for i in range(_n_dust)]
                elif len(self._minimum_specific_energy) != _n_dust:
                    raise Exception("Number of minimum_specific_energy values should match number of dust types")
                else:
                    _minimum_specific_energy = self._minimum_specific_energy

            else:

                _minimum_specific_energy = None

            if isinstance(self.grid, GridOnDisk):
                if _minimum_specific_energy is not None:
                    raise ValueError("Cannot set minimum specific energy or temperature when using grid from disk")
            elif _minimum_specific_energy is not None:
                g_grid['Quantities'].attrs["minimum_specific_energy"] = [float(x) for x in _minimum_specific_energy]

        else:

            root.create_group('Dust')

        root.close()

        self.filename = filename

    def add_point_source(self, *args, **kwargs):
        source = PointSource(*args, **kwargs)
        self.add_source(source)
        return source

    def add_point_source_collection(self, *args, **kwargs):
        source = PointSourceCollection(*args, **kwargs)
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

    def set_voronoi_grid(self, x, y, z,
                         xmin=None, xmax=None,
                         ymin=None, ymax=None,
                         zmin=None, zmax=None):
        self.set_grid(VoronoiGrid(x, y, z,
                                  xmin=xmin, xmax=xmax,
                                  ymin=ymin, ymax=ymax,
                                  zmin=zmin, zmax=zmax))

    def set_grid(self, grid):
        if isinstance(grid, AMRGrid):
            self.grid = AMRGrid(grid)
        else:
            if len(grid.quantities) > 0:
                self.grid = grid.__class__(grid)
            else:
                self.grid = deepcopy(grid)

    def use_grid_from_file(self, filename, path='/', dust=[]):
        """
        Use a grid from disk and don't read it in.

        Parameters
        ----------
        filename : str
            The name of the file to read from
        path : str, optional
            The path where the grid is located inside the file
        dust : list, optional
            A list containing a dust file or object for each dust index in the
            grid
        """
        self.grid = GridOnDisk(filename, path=path)
        self.dust = dust

    def add_peeled_images(self, sed=True, image=True):
        """
        Define a set of (peeled) images/SEDs.

        This returns a :class:`~hyperion.conf.PeeledImageConf` instance that
        can be used to set the image/SED parameters.

        Parameters
        ----------
        sed : bool, optional
            Whether to compute a set of SEDs
        image : bool, optional
            Whether to compute a set of images

        Returns
        -------
        image_conf : :class:`~hyperion.conf.PeeledImageConf` instance
            Instance that can be used to set the image/SED parameters

        Examples
        --------

        >>> image = m.add_peeled_images(sed=False, image=True)
        >>> image.set_wavelength_range(250, 0.01, 1000.)
        >>> ...
        """
        self.peeled_output.append(PeeledImageConf(sed=sed, image=image))
        self.peeled_output[-1]._set_monochromatic(self._monochromatic, frequencies=self._frequencies)
        return self.peeled_output[-1]

    def add_binned_images(self, **kwargs):
        """
        Define a set of (binned) images/SEDs.

        This returns a :class:`~hyperion.conf.BinnedImageConf` instance that
        can be used to set the image/SED parameters.

        Parameters
        ----------
        sed : bool, optional
            Whether to compute a set of SEDs
        image : bool, optional
            Whether to compute a set of images

        Returns
        -------
        image_conf : :class:`~hyperion.conf.BinnedImageConf` instance
            Instance that can be used to set the image/SED parameters

        Examples
        --------

        >>> image = m.add_binned_images(sed=False, image=True)
        >>> image.set_wavelength_range(250, 0.01, 1000.)
        >>> ...
        """
        if self.binned_output:
            raise Exception("Only one set of binned images can be set at this time")
        elif self._monochromatic:
            raise Exception("Binned images cannot be computed in monochromatic mode")
        else:
            self.binned_output = BinnedImageConf(sed=sed, image=image)
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
        """
        Run the model (should be called after `write()`).

        Parameters
        ----------
        filename : str, optional
            The output filename for the model. If not specified, then if the
            input file name contains ``.rtin``, then this is replaced with
            ``.rtout``, and otherwise ``.rtout`` is appended to the input
            filename.
        logfile : str, optional
            If specified, the standard output and errors will be output to
            this log file
        mpi : bool, optional
            Whether to run the model using the parallel (MPI) version of
            Hyperion.
        n_processes : int, optional
            If ``mpi`` is set to ``True``, this can be used to specify the
            number of processes to run Hyperion on.
        overwrite : bool, optional
            If set to ``True``, the output file is overwritten without
            warning.
        """

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
