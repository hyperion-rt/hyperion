import os
import string
import subprocess
import warnings

import h5py
import numpy as np
import matplotlib.pyplot as plt

import hyperion
from hyperion.util.functions import delete_file, random_id
from hyperion.grid import CartesianGrid, SphericalPolarGrid, CylindricalPolarGrid, OcTreeGrid, AMRGrid
from hyperion.sources import PointSource, SphericalSource, ExternalSphericalSource, ExternalBoxSource, MapSource, PlaneParallelSource
from hyperion.conf import RunConf, PeeledImageConf, BinnedImageConf, OutputConf
from hyperion.util.constants import c, pi
from hyperion.util.functions import FreezableClass
from hyperion.dust import SphericalDust

STOKESD = {}
STOKESD['I'] = 0
STOKESD['Q'] = 1
STOKESD['U'] = 2
STOKESD['V'] = 3

LABEL = {}
LABEL['I'] = '$\lambda\, F_\lambda$ (ergs/cm$^2$/s)'
LABEL['Q'] = '$\lambda\, F_\lambda$ [stokes=Q] (ergs/cm$^2$/s)'
LABEL['U'] = '$\lambda\, F_\lambda$ [stokes=U] (ergs/cm$^2$/s)'
LABEL['V'] = '$\lambda\, F_\lambda$ [stokes=V] (ergs/cm$^2$/s)'
LABEL['linpol'] = "Total linear polarization fraction"
LABEL['circpol'] = "Total circular polarization fraction"


def find_last_iteration(file_handle):
    max_iteration = 0
    for group_name in file_handle:
        if "Iteration" in group_name:
            iteration = int(group_name.split()[1])
            max_iteration = max(iteration, max_iteration)
    return max_iteration


def mc_linear_polarization(I, sigma_I, Q, sigma_Q, U, sigma_U):
    n = 1000
    if I > 0.:
        if sigma_I == 0.:
            Is = np.repeat(I, n)
        else:
            Is = np.random.normal(loc=I, scale=sigma_I, size=n)
        if sigma_Q == 0.:
            Qs = np.repeat(Q, n)
        else:
            Qs = np.random.normal(loc=Q, scale=sigma_Q, size=n)
        if sigma_U == 0.:
            Us = np.repeat(U, n)
        else:
            Us = np.random.normal(loc=U, scale=sigma_U, size=n)
        Ps = np.sqrt((Qs ** 2 + Us ** 2) / Is ** 2)
        return np.mean(Ps), np.std(Ps)
    else:
        return 0., 0.


def mc_circular_polarization(I, sigma_I, V, sigma_V):
    n = 1000
    if I > 0.:
        if sigma_I == 0.:
            Is = np.repeat(I, n)
        else:
            Is = np.random.normal(loc=I, scale=sigma_I, size=n)
        if sigma_V == 0.:
            Vs = np.repeat(V, n)
        else:
            Vs = np.random.normal(loc=V, scale=sigma_V, size=n)
        Ps = np.abs(Vs) / Is
        return np.mean(Ps), np.std(Ps)
    else:
        return 0., 0.


def bool2str(value):
    if value:
        return "yes"
    else:
        return "no"


class File(file):

    def __init__(self, filename, *args, **kwargs):
        self.abspath = os.path.abspath(filename)
        file.__init__(self, filename, *args, **kwargs)


class Configuration(FreezableClass):

    def __init__(self):
        self.output = OutputConf()
        self.run = RunConf()
        self._freeze()


class Model(FreezableClass):

    def __init__(self, name):
        '''
        Create a Model instance

        Parameters
        ----------
        name : str
            The name of the model. This is used to name the .rtin file
            output with the write() method.
        '''

        if name[-6:] == '.rtout':
            name = name[:-6]

        self.conf = Configuration()

        self.name = name
        self.reset_density()
        self.reset_sources()
        self.reset_images()
        self.specific_energy = None
        self.minimum_specific_energy = []
        self.set_monochromatic(False)

        self.grid = None

        # Import methods for convenience
        self.set_n_initial_iterations = self.conf.run.set_n_initial_iterations
        self.set_n_photons = self.conf.run.set_n_photons
        self.set_raytracing = self.conf.run.set_raytracing
        self.set_max_interactions = self.conf.run.set_max_interactions
        self.set_pda = self.conf.run.set_pda
        self.set_mrw = self.conf.run.set_mrw
        self.set_convergence = self.conf.run.set_convergence
        self.set_kill_on_absorb = self.conf.run.set_kill_on_absorb
        self.set_forced_first_scattering = self.conf.run.set_forced_first_scattering
        self.set_output_bytes = self.conf.run.set_output_bytes
        self.set_sample_sources_evenly = self.conf.run.set_sample_sources_evenly
        self.set_enforce_energy_range = self.conf.run.set_enforce_energy_range

        self._freeze()

    def reset_density(self):
        self.density = []
        self.dust = []

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
        self.conf.run._monochromatic = monochromatic

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
            group.create_dataset('Frequencies', data=np.array(zip(self._frequencies), dtype=[('nu', dtype)]), compression=compression)

    def write(self, compression=True, copy_dust=False, absolute_paths=False, geo_dtype=float, wall_dtype=float, physics_dtype=float):
        '''
        Write the model to an HDF5 file

        Optional Keyword Arguments:

            *compression*: [ True | False ]
                Whether to compress the datasets inside the HDF5 file

            *copy_dust*: [ True | False ]
                Whether to copy the dust properties when a filename is specified (True), or to just link to the file (False)

            *paths*: [ True | False ]
                Whether to use absolute paths (True) or relative paths (False) when linking to files
        '''

        if not self.name:
            raise Exception("Model name is not defined")

        if not os.path.dirname(self.name) == "" and not os.path.exists(os.path.dirname(self.name)):
            raise Exception("Directory %s does not exist" % os.path.dirname(self.name))

        # Create output file
        delete_file('%s.rtin' % self.name)
        root = h5py.File('%s.rtin' % self.name, 'w')
        root.attrs['python_version'] = hyperion.__version__

        g_dust = root.create_group('Dust')
        g_grid = root.create_group('Grid')
        g_geometry = g_grid.create_group('Geometry')
        g_physics = g_grid.create_group('Physics')
        g_sources = root.create_group('Sources')
        g_output = root.create_group('Output')
        g_peeled = g_output.create_group('Peeled')
        g_binned = g_output.create_group('Binned')

        if self.grid is None:
            raise Exception("No grid has been set up")

        # Generate random geometry ID
        self.grid.geometry_id = random_id()

        # Output sources
        for i, source in enumerate(self.sources):
            if isinstance(source, MapSource):
                source.write(g_sources, 'Source %05i' % i, self.grid, compression=compression, map_dtype=physics_dtype)
            else:
                source.write(g_sources, 'Source %05i' % i)

        # Output configuration files
        for i, peel in enumerate(self.peeled_output):
            if not self._frequencies is None:
                if not peel._monochromatic:
                    raise Exception("Peeled images need to be set to monochromatic mode")
            peel.write(g_peeled.create_group('Group %05i' % (i + 1)))
        if self.binned_output is not None:
            self.binned_output.write(g_binned.create_group('Group 00001'))

        self._write_monochromatic(root, compression=compression)
        self.conf.run.write(root)
        self.conf.output.write(g_output)

        if len(self.density) > 0:

            self.grid.write_physical_array(g_physics, self.density, "Density", dust=True, compression=compression, physics_dtype=physics_dtype)
            g_physics.create_dataset("Minimum Specific Energy", data=self.minimum_specific_energy)

            if self.specific_energy is not None:
                if type(self.specific_energy) is list:
                    self.grid.write_physical_array(g_physics, self.specific_energy, "Specific Energy", dust=True, compression=compression, physics_dtype=physics_dtype)
                elif type(self.specific_energy) is str:
                    f = h5py.File(self.specific_energy, 'r')
                    max_iteration = find_last_iteration(f)
                    print "Retrieving specific_energy from iteration %i of %s" % (max_iteration, self.specific_energy)
                    specific_energy = f['Iteration %05i' % max_iteration]['specific_energy']
                    if specific_energy.shape != g_physics['Density'].shape:
                        raise Exception("Specific energy array is not the right dimensions")
                    dset = g_physics.create_dataset('Specific Energy', data=specific_energy)
                    dset.attrs['geometry'] = g_physics['Density'].attrs['geometry']
                else:
                    raise Exception("Unknown type %s for Model.specific_energy" % str(type(self.specific_energy)))

            # Output dust file
            present = {}
            for i, dust in enumerate(self.dust):

                short_name = 'dust_%03i' % (i + 1)

                if copy_dust:

                    if type(dust) == str:
                        dust = SphericalDust(dust)

                    if id(dust) in present:
                        g_dust[short_name] = h5py.SoftLink(present[id(dust)])
                    else:
                        dust.write(g_dust.create_group(short_name))
                        present[id(dust)] = short_name

                else:

                    if type(dust) != str:
                        dust = dust.filename

                    if absolute_paths:
                        g_dust[short_name] = h5py.ExternalLink(os.path.abspath(dust), '/')
                    else:
                        g_dust[short_name] = h5py.ExternalLink(os.path.relpath(dust), '/')

        # Output geometry
        self.grid.write_geometry(g_geometry, \
                                 compression=compression, geo_dtype=geo_dtype, wall_dtype=wall_dtype)

        root.close()

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

    def add_density_grid(self, density, dust, specific_energy=None, minimum_specific_energy=None, minimum_temperature=None, merge_if_possible=True):

        # Check that grid has been previously defined
        if not self.grid:
            raise Exception("Grid not defined")

        # Check whether grid dimensions are the same
        if not isinstance(self.grid, AMRGrid):
            if not density.shape == self.grid.shape:
                raise Exception("Density shape does not match that of grid")
            if specific_energy is not None:
                if not specific_energy.shape == self.grid.shape:
                    raise Exception("Specific Energy shape does not match that of grid")

        # Check whether all densities are zero
        if np.all(density == 0.):
            warnings.warn("All density values are zero - ignoring density grid")
            return

        # Check consistency between density list size and specific energy list size
        if len(self.density) > 0:
            if specific_energy is not None and type(self.specific_energy) is not list:
                raise Exception("Cannot add specific energy as it was not added for previous density arrays")
            if specific_energy is None and type(self.specific_energy) is list:
                raise Exception("Specific energy was added for previous density arrays, so should be added for all arrays")
        else:
            if specific_energy is not None:
                self.specific_energy = []

        if minimum_specific_energy is not None and minimum_temperature is not None:
            raise Exception("Cannot specify both the minimum specific energy and temperature")
        elif minimum_temperature is not None:
            d = SphericalDust(dust) if type(dust) == str else dust
            minimum_specific_energy = d.mean_opacities._temperature2specific_energy(minimum_temperature)

        # Check whether the density can be added to an existing one
        if merge_if_possible:

            # Only consider this if the specific energy is not specified
            if specific_energy is None:

                # Only do it if the dust type already exists
                if dust in self.dust:

                    ip = self.dust.index(dust)

                    # Check whether the minimum_specific_energy values differ
                    if minimum_specific_energy is not None:
                        if self.minimum_specific_energy[ip] != minimum_specific_energy:
                            warnings.warn("Cannot merge density grids because minimum_specific_energy values differ")
                            merge = False
                        else:
                            merge = True
                    else:
                        merge = True

                    # Merge the densities
                    if merge:
                        warnings.warn("Merging densities")
                        self.density[ip] += density
                        return

        # Set the density and dust
        self.density.append(density)
        self.dust.append(dust)

        # Set specific energy if specified
        if specific_energy is not None:
            self.specific_energy.append(specific_energy)

        # Set minimum specific energy
        if minimum_specific_energy is not None:
            self.minimum_specific_energy.append(minimum_specific_energy)
        else:
            self.minimum_specific_energy.append(0.)

    def set_cartesian_grid(self, x_wall, y_wall, z_wall):
        self.grid = CartesianGrid(x_wall, y_wall, z_wall)

    def set_cylindrical_polar_grid(self, w_wall, z_wall, p_wall):
        self.grid = CylindricalPolarGrid(w_wall, z_wall, p_wall)

    def set_spherical_polar_grid(self, r_wall, t_wall, p_wall):
        self.grid = SphericalPolarGrid(r_wall, t_wall, p_wall)

    def set_octree_grid(self, refined, x, y, z, dx, dy, dz):
        self.grid = OcTreeGrid(refined, x, y, z, dx, dy, dz)

    def set_amr_grid(self, description):
        self.grid = AMRGrid(description)

    def add_peeled_images(self, **kwargs):
        self.peeled_output.append(PeeledImageConf(**kwargs))
        return self.peeled_output[-1]

    def add_binned_images(self, **kwargs):
        if self.binned_output:
            raise Exception("Only one set of binned images can be set at this time")
        else:
            self.binned_output = BinnedImageConf(**kwargs)
            return self.binned_output

    def run(self, logfile=None):
        input = "%s.rtin" % self.name
        output = "%s.rtout" % self.name
        if logfile:
            flog = file(logfile, 'wb')
            p = subprocess.call('hyperion %s %s' % (input, output), stdout=flog, stderr=flog, shell=True)
        else:
            p = subprocess.call('hyperion %s %s' % (input, output), shell=True)

    def get_sed(self, stokes='I', group=1, technique='peeled',
                distance=None, component='total', inclination='all',
                aperture='all', uncertainties=False, units='ergs/cm^2/s',
                source_id=None, dust_id=None):
        '''
        Retrieve SEDs for a specific image group and Stokes component

        Parameters
        ----------

        stokes : str, optional
            The Stokes component to return. This can be:
                * 'I': Total intensity [default]
                * 'Q': Q Stokes parameter (linear polarization)
                * 'U': U Stokes parameter (linear polarization)
                * 'V': V Stokes parameter (circular polarization)
                * 'linpol':  Total linear polarization fraction
                * 'circpol': Total circular polariation fraction

        technique : str, optional
            Whether to retrieve SED(s) computed with photon peeling-off
            ('peeled') or binning ('binned'). Default is 'peeled'.

        group : int, optional
            The peeloff group. If multiple peeloff image groups were required,
            this can be used to select between them. The default is to return
            the first group. This option is only used if technique='peeled'.

        distance : float, optional
            The distance to the observer, in cm.

        component : str, optional
            The component to return based on origin and last interaction.
            This can be:

                * 'total': Total flux

                * 'source_emit': The photons were last emitted from a source
                  and did not undergo any subsequent interactions.

                * 'dust_emit': The photons were last emitted dust and did not
                  undergo any subsequent interactions

                * 'source_scat': The photons were last emitted from a source
                  and were subsequently scattered

                * 'dust_scat': The photons were last emitted from dust and
                  were subsequently scattered

        aperture : int, optional
            The number of the aperture to plot (zero-based). Use 'all' to
            return all apertures, and -1 to show the largest aperture.

        inclination : int, optional
            The number of the viewing angle to plot (zero-based). Use 'all'
            to return all viewing angles.

        uncertainties : bool, optional
            Whether to compute and return uncertainties

        units : str, optional
            The output units for the SED(s). Valid options if a distance is
            specified are:
                * ``'ergs/cm^2/s'``
                * ``'ergs/cm^2/s/Hz'``
                * ``'Jy'``
                * ``'mJy'``
                * ``'MJy/sr'``
            The default is ``'ergs/cm^2/s'``. If a distance is not specified,
            then this option is ignored, and the output units are ergs/s.

        source_id, dust_id : int, optional
            If the output file was made with track_origin='detailed', a
            specific source and dust component can be specified. If
            these are not specified, then the total component requested
            for all sources or dust types is returned.

        Returns
        -------

        wav : numpy.ndarray
            The wavelengths for which the SEDs are defined, in microns

        flux or degree of polarization : numpy.ndarray
            The flux or degree of polarization. This is a data cube which has
            at most three dimensions (n_inclinations, n_apertures,
            n_wavelengths). If an aperture or inclination is specified, this
            reduces the number of dimensions in the flux cube. If `stokes` is
            one of 'I', 'Q', 'U', or 'V', the flux is either returned in
            ergs/s (if distance is not specified) or in the units specified by
            units= (if distance is specified). If `stokes` is one of 'linpol'
            or 'circpol', the degree of polarization is returned as a fraction
            in the range 0 to 1.

        uncertainty : numpy.ndarray (if `uncertainties`=True)
            The uncertainties on the flux or degree of polarization. This
            has the same dimensions as the flux or degree of polarization
            array.
        '''

        # Check for inconsistent parameters
        if distance is not None and stokes in ['linpol', 'circpol']:
            raise Exception("Cannot scale linear or circular polarization degree by distance")

        # Check that file exists
        if not os.path.exists('%s.rtout' % self.name):
            raise Exception("File not found %s.rtout" % self.name)

        # Open file and retrieve correct group
        f = h5py.File('%s.rtout' % self.name, 'r')

        if technique == 'peeled':
            g = f['Peeled/Group %05i' % group]
        else:
            g = f['Binned']

        # Check that uncertainties are present if requested
        if uncertainties and not 'seds_unc' in g:
            raise Exception("Uncertainties requested but not present in file")

        if 'track_origin' in g['seds'].attrs:

            track_origin = g['seds'].attrs['track_origin']

            if track_origin == 'no' and component != 'total':
                raise Exception("cannot extract component=%s - file only contains total flux" % component)

            if track_origin != 'detailed':
                if source_id is not None:
                    raise Exception("cannot specify source_id, as SEDs were not computed with track_origin='detailed'")
                if dust_id is not None:
                    raise Exception("cannot specify dust_id, as SEDs were not computed with track_origin='detailed'")

            if component in ['source_emit', 'dust_emit', 'source_scat', 'dust_scat']:

                if component == 'source_emit':
                    io = 0
                elif component == 'dust_emit':
                    io = 1
                elif component == 'source_scat':
                    io = 2
                elif component == 'dust_scat':
                    io = 3

                if track_origin == 'detailed':

                    ns = g['seds'].attrs['n_sources']
                    nd = g['seds'].attrs['n_dust']

                    io = ((io - io % 2) * ns + (io - (io + 1) % 2 - 1) * nd) / 2

                    if component.startswith('source'):
                        if source_id is None:
                            io = (io, io + ns)
                        else:
                            if source_id < 1 or source_id > ns:
                                raise Exception("source_id should be between 1 and %i" % ns)
                            io = io + (source_id - 1)
                    else:
                        if dust_id is None:
                            io = (io, io + nd)
                        else:
                            if dust_id < 1 or dust_id > nd:
                                raise Exception("dust_id should be between 1 and %i" % nd)
                            io = io + (dust_id - 1)

        # Set up wavelength space
        if 'numin' in g['seds'].attrs:
            numin = g['seds'].attrs['numin']
            numax = g['seds'].attrs['numax']
            wavmin, wavmax = c / numax * 1.e4, c / numin * 1.e4
            wav = np.logspace(np.log10(wavmax), np.log10(wavmin), g['seds'].shape[4] * 2 + 1)[1::2]
            nu = c / wav * 1.e4
        else:
            nu = g['frequencies']['nu']
            wav = c / nu * 1.e4

        # Optionally scale by distance
        if distance:
            scale = 1. / (4. * pi * distance ** 2)
        else:
            scale = 1.

        seds = g['seds'].value
        if uncertainties:
            seds_unc = g['seds_unc'].value

        try:
            inside_observer = g.attrs['inside_observer']
        except:
            inside_observer = False

        # Optionally scale by distance
        if distance or inside_observer:

            # Convert to the correct units
            if units == 'ergs/cm^2/s':
                scale = np.repeat(1., len(nu))
            elif units == 'ergs/cm^2/s/Hz':
                scale = 1. / nu
            elif units == 'Jy':
                scale = 1.e23 / nu
            elif units == 'mJy':
                scale = 1.e26 / nu
            elif units == 'MJy/sr':

                # Find spatial extent of the image
                xmin = g['images'].attrs['xmin']
                xmax = g['images'].attrs['xmax']
                ymin = g['images'].attrs['ymin']
                ymax = g['images'].attrs['ymax']

                # Find pixel dimensions of image
                ny, nx = images.shape[-2:]

                # Find pixel resolution in radians/pixel
                if inside_observer:
                    pix_dx = abs(np.radians(xmax - xmin)  / float(nx))
                    pix_dy = abs(np.radians(ymax - ymin)  / float(ny))
                else:
                    pix_dx = abs(np.arctan((xmax - xmin) / float(nx) / distance))
                    pix_dy = abs(np.arctan((ymax - ymin) / float(ny) / distance))

                # Find pixel area in steradians
                pix_area = pix_dx * pix_dy

                scale = 1.e17 / nu / pix_area

            else:
                raise Exception("Unknown units: %s" % scale)

            # Scale by distance
            if distance:
                scale *= 1. / (4. * pi * distance ** 2)

        else:

            # Units here are not technically ergs/cm^2/s but ergs/s
            scale = np.repeat(1., len(nu))

        # If in 32-bit mode, need to convert to 64-bit because of scaling/polarization to be safe
        if seds.dtype == np.float32:
            seds = seds.astype(np.float64)
        if uncertainties and seds_unc.dtype == np.float32:
            seds_unc = seds_unc.astype(np.float64)

        # Select correct origin component
        if component == 'total':
            flux = np.sum(seds[:, :, :, :, :], axis=1)
            if uncertainties:
                unc = np.sqrt(np.sum(seds_unc[:, :, :, :, :] ** 2, axis=1))
        elif component in ['source_emit', 'dust_emit', 'source_scat', 'dust_scat']:
            if type(io) is tuple:
                start, end = io
                flux = np.sum(seds[:, start:end, :, :, :], axis=1)
                if uncertainties:
                    unc = np.sqrt(np.sum(seds_unc[:, start:end, :, :, :] ** 2, axis=1))
            else:
                flux = seds[:, io, :, :, :]
                if uncertainties:
                    unc = seds_unc[:, io, :, :, :]
        else:
            raise Exception("Unknown component: %s" % component)

        # Close HDF5 file
        f.close()

        # Select correct Stokes component
        if stokes in STOKESD:
            y = flux[STOKESD[stokes], :, :, :] * scale[np.newaxis, np.newaxis, :]
            if uncertainties:
                yerr = unc[STOKESD[stokes], :, :, :] * scale[np.newaxis, np.newaxis, :]
            else:
                yerr = np.zeros(y.shape)
        elif stokes == 'linpol':
            if uncertainties:
                f = np.vectorize(mc_linear_polarization)
                y, yerr = f(flux[0, :, :, :], unc[0, :, :, :], flux[1, :, :, :], unc[1, :, :, :], flux[2, :, :, :], unc[2, :, :, :])
            else:
                y = np.sqrt((flux[1, :, :, :] ** 2 + flux[2, :, :, :] ** 2) / flux[0, :, :, :] ** 2)
                y[np.isnan(y)] = 0.
                yerr = np.zeros(y.shape)
        elif stokes == 'circpol':
            if uncertainties:
                f = np.vectorize(mc_circular_polarization)
                y, yerr = f(flux[0, :, :, :], unc[0, :, :, :], flux[3, :, :, :], unc[3, :, :, :])
            else:
                y = np.abs(flux[3, :, :, :] / flux[0, :, :, :])
                y[np.isnan(y)] = 0.
                yerr = np.zeros(y.shape)
        else:
            raise Exception("Unknown Stokes parameter: %s" % stokes)

        if inclination == 'all' and aperture == 'all':
            pass
        elif inclination == 'all':
            y, yerr = y[:, aperture, :], yerr[:, aperture, :]
        elif aperture == 'all':
            y, yerr = y[inclination, :, :], yerr[inclination, :, :]
        else:
            y, yerr = y[inclination, :, :], yerr[:, aperture, :]
            if np.isscalar(inclination):
                y, yerr = y[aperture, :], yerr[aperture, :]
            else:
                y, yerr = y[:, aperture, :], yerr[:, aperture, :]

        if uncertainties:
            return wav, y, yerr
        else:
            return wav, y

    def plot_sed(self, axes=None, filename=None, stokes='I',
                     technique='peeled', group=1, distance=None,
                     component='total', aperture=-1, inclination='all',
                     wmin=0.01, wmax=5000., fmin=None, fmax=None,
                     color='black', labels=True, uncertainties=False,
                     source_id=None, dust_id=None, **kwargs):
        '''
        Plot an SED or polarization spectrum

        Parameters
        ----------

        axes : matplotlib.pyplot.Axes instance, optional
            The matplotlib Axes to plot the SED(s) into (incompatible with
            the `filename` option).

        filename : str, optional
            The file to write the plot to (incompatible with the `ax`
            option).

        stokes : str, optional
            The Stokes component to plot. This can be:
                'I'       Total intensity [default]
                'Q'       Q Stokes parameter (linear polarization)
                'U'       U Stokes parameter (linear polarization)
                'V'       V Stokes parameter (circular polarization)
                'linpol'  Total linear polarization fraction
                'circpol' Total circular polariation fraction

        technique : str, optional
            Whether to plot SED(s) computed with photon peeling-off
            ('peeled') or binning ('binned'). Default is 'peeled'.

        group : int, optional
            The peeloff group. If multiple peeloff image groups were required,
            this can be used to select between them. The default is to return
            the first group. This option is only used if technique='peeled'.

        distance : float, optional
            The distance to the observer, in cm.

        component : str, optional
            The component to plot based on origin and last interaction.
            This can be:
                'total'       Total SED
                'source_emit' The photons were last emitted from a
                              source and did not undergo any subsequent
                              interactions.
                'dust_emit'   The photons were last emitted dust and did
                              not undergo any subsequent interactions
                'source_scat' The photons were last emitted from a
                              source and were subsequently scattered
                'dust_scat'   The photons were last emitted from dust and
                              were subsequently scattered

        aperture : int
            The number of the aperture to plot (zero-based). Use 'all' to
            show all apertures, and -1 to show the largest aperture.

        inclination : int, optional
            The number of the viewing angle to plot (zero-based). Use 'all'
            to show all viewing angles.

        wmin, wmax : float, optional
            The range in wavelengths to show (in microns).

        fmin, fmax : float, optional
            The range in fluxes to show (in ergs/s if distance is not
            specified, or in ergs/cm^2/s if distance is specified).

        color : str, optional
            The color to plot the SED(s) in.

        labels : bool, optional
            Whether or not to show the axis labels

        uncertainties : bool, optional
            Whether or not to show the uncertainties

        source_id, dust_id : int, optional
            If the output file was made with track_origin='detailed', a
            specific source and dust component can be specified. If
            these are not specified, then the total component requested
            for all sources or dust types is returned.

        Returns
        -------

        axes : matplotlib.pyplot.Axes instance (if `axes` was set)
            The updated matplotlib Axes

        figure : matplotlib.pyplot.Figure instance (if `axes` and `filename` were not set)
            The matplotlib Figure created

        '''

        # Check for inconsistent parameters
        if axes is not None and filename is not None:
            raise Exception("Cannot specify both an axes instance and a filename")

        if axes is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = axes

        # Make sure that SED cube always has three dimensions, so need to request inclinations and apertures in a list
        if type(inclination) is not str and np.isscalar(inclination):
            inclination = [inclination]
        if type(aperture) is not str and np.isscalar(aperture):
            aperture = [aperture]

        seds = self.get_sed(stokes=stokes, technique=technique, group=group, distance=distance, component=component, aperture=aperture, inclination=inclination, uncertainties=uncertainties, source_id=source_id, dust_id=dust_id)

        if uncertainties:
            wav, nufnu, nufnu_unc = seds
        else:
            wav, nufnu = seds

        if inclination == 'all':
            inclination = range(nufnu.shape[0])

        if aperture == 'all':
            aperture = range(nufnu.shape[1])

        maxflux = 0.
        for ii in inclination:
            for ia in aperture:
                if(np.all(nufnu[ii, ia, :] == 0.)):
                    warnings.warn("All fluxes are zero - cannot plot in log-log plot")
                else:

                    # Plot main SED
                    if stokes in ['Q', 'U', 'V']:
                        ax.plot(wav, nufnu[ii, ia, :], color=color, **kwargs)
                        ax.set_xscale('log')
                    else:
                        ax.loglog(wav, nufnu[ii, ia, :], color=color, **kwargs)

                    # Check whether the maximum flux needs to be changed
                    maxflux = max(maxflux, np.nanmax(nufnu))

                    if uncertainties:
                        nufnu_lower = np.maximum(nufnu[ii, ia, :] * 1.e-10, nufnu[ii, ia, :] - nufnu_unc[ii, ia, :])
                        nufnu_unc_lower = nufnu[ii, ia, :] - nufnu_lower
                        ax.errorbar(wav, nufnu[ii, ia, :], yerr=[nufnu_unc_lower, nufnu_unc[ii, ia, :]], fmt=None, ecolor=color, **kwargs)

        if not fmax:
            if stokes in ['Q', 'U', 'V']:
                fmax = maxflux * 1.5
            elif stokes in ['linpol', 'circpol']:
                fmax = 1.
            else:
                fmax = maxflux * 10.

        if not fmin:
            if stokes in ['Q', 'U', 'V']:
                fmin = -fmax
            else:
                fmin = fmax * 1.e-6

        ax.set_xlim(wmin, wmax)
        ax.set_ylim(fmin, fmax)

        if labels:
            ax.set_xlabel('$\lambda$ ($\mu$m)')
            ax.set_ylabel(LABEL[stokes])

        if filename:
            fig.savefig(filename)
            plt.close(fig)
            return
        elif axes is None:
            return fig
        else:
            return ax

    def get_image(self, stokes='I', group=1, technique='peeled',
                  distance=None, component='total', inclination='all',
                  uncertainties=False, units='ergs/cm^2/s',
                  source_id=None, dust_id=None):
        '''
        Retrieve images for a specific image group and Stokes component

        Parameters
        ----------

        stokes : str, optional
            The Stokes component to return. This can be:
                * 'I': Total intensity [default]
                * 'Q': Q Stokes parameter (linear polarization)
                * 'U': U Stokes parameter (linear polarization)
                * 'V': V Stokes parameter (circular polarization)
                * 'linpol':  Total linear polarization fraction
                * 'circpol': Total circular polariation fraction

        technique : str, optional
            Whether to retrieve an image computed with photon peeling-off
            ('peeled') or binning ('binned'). Default is 'peeled'.

        group : int, optional
            The peeloff group. If multiple peeloff image groups were required,
            this can be used to select between them. The default is to return
            the first group. This option is only used if technique='peeled'.

        distance : float, optional
            The distance to the observer, in cm.

        component : str, optional
            The component to return based on origin and last interaction.
            This can be:

                * 'total': Total flux

                * 'source_emit': The photons were last emitted from a source
                  and did not undergo any subsequent interactions.

                * 'dust_emit': The photons were last emitted dust and did not
                  undergo any subsequent interactions

                * 'source_scat': The photons were last emitted from a source
                  and were subsequently scattered

                * 'dust_scat': The photons were last emitted from dust and
                  were subsequently scattered

        inclination : int, optional
            The number of the viewing angle to plot (zero-based). Use 'all'
            to return all viewing angles.

        uncertainties : bool, optional
            Whether to compute and return uncertainties

        units : str, optional
            The output units for the image(s). Valid options if a distance is
            specified are:
                * ``'ergs/cm^2/s'``
                * ``'ergs/cm^2/s/Hz'``
                * ``'Jy'``
                * ``'mJy'``
                * ``'MJy/sr'``
            The default is ``'ergs/cm^2/s'``. If a distance is not specified,
            then this option is ignored, and the output units are ergs/s.

        source_id, dust_id : int, optional
            If the output file was made with track_origin='detailed', a
            specific source and dust component can be specified. If
            these are not specified, then the total component requested
            for all sources or dust types is returned.

        Returns
        -------

        wav : numpy.ndarray
            The wavelengths for which the SEDs are defined, in microns

        flux or degree of polarization : numpy.ndarray
            The flux or degree of polarization. This is a data cube which has
            at most three dimensions (n_inclinations, n_wavelengths). If an
            aperture or inclination is specified, this reduces the number of
            dimensions in the flux cube. If `stokes` is one of 'I', 'Q', 'U',
            or 'V', the flux is either returned in ergs/s (if distance is not
            specified) or in the units specified by units= (if distance is
            specified). If `stokes` is one of 'linpol' or 'circpol', the
            degree of polarization is returned as a fraction in the range 0 to
            1.

        uncertainty : numpy.ndarray (if `uncertainties`=True)
            The uncertainties on the flux or degree of polarization. This
            has the same dimensions as the flux or degree of polarization
            array.
        '''

        # Check for inconsistent parameters
        if distance is not None and stokes in ['linpol', 'circpol']:
            raise Exception("Cannot scale linear or circular polarization degree by distance")

        # Check that file exists
        if not os.path.exists('%s.rtout' % self.name):
            raise Exception("File not found %s.rtout" % self.name)

        # Open file and retrieve correct group
        f = h5py.File('%s.rtout' % self.name, 'r')

        if technique == 'peeled':
            g = f['Peeled/Group %05i' % group]
        else:
            g = f['Binned']

        # Check that uncertainties are present if requested
        if uncertainties and not 'images_unc' in g:
            raise Exception("Uncertainties requested but not present in file")

        if 'track_origin' in g['images'].attrs:

            track_origin = g['images'].attrs['track_origin']

            if track_origin == 'no' and component != 'total':
                raise Exception("cannot extract component=%s - file only contains total flux" % component)

            if track_origin != 'detailed':
                if source_id is not None:
                    raise Exception("cannot specify source_id, as images were not computed with track_origin='detailed'")
                if dust_id is not None:
                    raise Exception("cannot specify dust_id, as images were not computed with track_origin='detailed'")

            if component in ['source_emit', 'dust_emit', 'source_scat', 'dust_scat']:

                if component == 'source_emit':
                    io = 0
                elif component == 'dust_emit':
                    io = 1
                elif component == 'source_scat':
                    io = 2
                elif component == 'dust_scat':
                    io = 3

                if track_origin == 'detailed':

                    ns = g['images'].attrs['n_sources']
                    nd = g['images'].attrs['n_dust']

                    io = ((io - (io + 1) % 2 + 1) * ns + (io - io % 2) * nd) / 2

                    if component.startswith('source'):
                        if source_id is None:
                            io = (io, io + ns)
                        else:
                            if source_id < 1 or source_id > ns:
                                raise Exception("source_id should be between 1 and %i" % ns)
                            io = io + (source_id - 1)
                    else:
                        if dust_id is None:
                            io = (io, io + nd)
                        else:
                            if dust_id < 1 or dust_id > nd:
                                raise Exception("dust_id should be between 1 and %i" % nd)
                            io = io + (dust_id - 1)

        # Set up wavelength space
        if 'numin' in g['images'].attrs:
            numin = g['images'].attrs['numin']
            numax = g['images'].attrs['numax']
            wavmin, wavmax = c / numax * 1.e4, c / numin * 1.e4
            wav = np.logspace(np.log10(wavmax), np.log10(wavmin), g['images'].shape[3] * 2 + 1)[1::2]
            nu = c / wav * 1.e4
        else:
            nu = g['frequencies']['nu']
            wav = c / nu * 1.e4

        images = g['images'].value
        if uncertainties:
            images_unc = g['images_unc'].value

        try:
            inside_observer = g.attrs['inside_observer']
        except:
            inside_observer = False

        # Optionally scale by distance
        if distance or inside_observer:

            # Convert to the correct units
            if units == 'ergs/cm^2/s':
                scale = np.repeat(1., len(nu))
            elif units == 'ergs/cm^2/s/Hz':
                scale = 1. / nu
            elif units == 'Jy':
                scale = 1.e23 / nu
            elif units == 'mJy':
                scale = 1.e26 / nu
            elif units == 'MJy/sr':

                # Find spatial extent of the image
                xmin = g['images'].attrs['xmin']
                xmax = g['images'].attrs['xmax']
                ymin = g['images'].attrs['ymin']
                ymax = g['images'].attrs['ymax']

                # Find pixel dimensions of image
                ny, nx = images.shape[-2:]

                # Find pixel resolution in radians/pixel
                if inside_observer:
                    pix_dx = abs(np.radians(xmax - xmin)  / float(nx))
                    pix_dy = abs(np.radians(ymax - ymin)  / float(ny))
                else:
                    pix_dx = abs(np.arctan((xmax - xmin) / float(nx) / distance))
                    pix_dy = abs(np.arctan((ymax - ymin) / float(ny) / distance))

                # Find pixel area in steradians
                pix_area = pix_dx * pix_dy

                scale = 1.e17 / nu / pix_area

            else:
                raise Exception("Unknown units: %s" % scale)

            # Scale by distance
            if distance:
                scale *= 1. / (4. * pi * distance ** 2)

        else:

            # Units here are not technically ergs/cm^2/s but ergs/s
            scale = np.repeat(1., len(nu))

        # If in 32-bit mode, need to convert to 64-bit because of scaling/polarization to be safe
        if images.dtype == np.float32:
            images = images.astype(np.float64)
        if uncertainties and images_unc.dtype == np.float32:
            images_unc = images_unc.astype(np.float64)

        # Select correct origin component
        if component == 'total':
            flux = np.sum(images[:, :, :, :, :, :], axis=1)
            if uncertainties:
                unc = np.sqrt(np.sum(images_unc[:, :, :, :, :, :] ** 2, axis=1))
        elif component in ['source_emit', 'dust_emit', 'source_scat', 'dust_scat']:
            if type(io) is tuple:
                start, end = io
                flux = np.sum(images[:, start:end, :, :, :, :], axis=1)
                if uncertainties:
                    unc = np.sqrt(np.sum(images_unc[:, start:end, :, :, :, :] ** 2, axis=1))
            else:
                flux = images[:, io, :, :, :, :]
                if uncertainties:
                    unc = images_unc[:, io, :, :, :, :]
        else:
            raise Exception("Unknown component: %s" % component)

        # Close HDF5 file
        f.close()

        # Select correct Stokes component
        if stokes in STOKESD:
            y = flux[STOKESD[stokes], :, :, :, :] * scale[np.newaxis, :, np.newaxis, np.newaxis]
            if uncertainties:
                yerr = unc[STOKESD[stokes], :, :, :, :] * scale[np.newaxis, :, np.newaxis, np.newaxis]
            else:
                yerr = np.zeros(y.shape)
        elif stokes == 'linpol':
            if uncertainties:
                f = np.vectorize(mc_linear_polarization)
                y, yerr = f(flux[0, :, :, :, :], unc[0, :, :, :, :], flux[1, :, :, :, :], unc[1, :, :, :, :], flux[2, :, :, :, :], unc[2, :, :, :, :])
            else:
                y = np.sqrt((flux[1, :, :, :, :] ** 2 + flux[2, :, :, :, :] ** 2) / flux[0, :, :, :, :] ** 2)
                y[np.isnan(y)] = 0.
                yerr = np.zeros(y.shape)
        elif stokes == 'circpol':
            if uncertainties:
                f = np.vectorize(mc_circular_polarization)
                y, yerr = f(flux[0, :, :, :, :], unc[0, :, :, :, :], flux[3, :, :, :, :], unc[3, :, :, :, :])
            else:
                y = np.abs(flux[3, :, :, :, :] / flux[0, :, :, :, :])
                y[np.isnan(y)] = 0.
                yerr = np.zeros(y.shape)
        else:
            raise Exception("Unknown Stokes parameter: %s" % stokes)

        if inclination == 'all':
            pass
        else:
            y, yerr = y[inclination, :, :, :], yerr[inclination, :, :, :]

        if uncertainties:
            return wav, y, yerr
        else:
            return wav, y

    def plot_image(self, wavelength, axes=None, filename=None, stokes='I', group=1,
                   technique='peeled', distance=None, component='total',
                   inclination=0, vmin=None, vmax=None, cmap=None,
                   labels=True, source_id=None, dust_id=None, **kwargs):
        '''
        Plot an image

        Parameters
        ----------

        axes : matplotlib.pyplot.Axes instance, optional
            The matplotlib Axes to plot the SED(s) into (incompatible with
            the `filename` option).

        wavelength : float
            The wavelength to plot the image for, in microns. The
            wavelength closest to this will be chose.

        stokes : str, optional
            The Stokes component to return. This can be:
                'I'       Total intensity [default]
                'Q'       Q Stokes parameter (linear polarization)
                'U'       U Stokes parameter (linear polarization)
                'V'       V Stokes parameter (circular polarization)
                'linpol'  Total linear polarization fraction
                'circpol' Total circular polariation fraction

        technique : str, optional
            Whether to plot an image computed with photon peeling-off
            ('peeled') or binning ('binned'). Default is 'peeled'.

        group : int, optional
            The peeloff group. If multiple peeloff image groups were required,
            this can be used to select between them. The default is to return
            the first group. This option is only used if technique='peeled'.

        distance : float, optional
            The distance to the observer, in cm.

        component : str, optional
            The component to return based on origin and last interaction.
            This can be:
                'total'       Total SED
                'source_emit' The photons were last emitted from a
                              source and did not undergo any subsequent
                              interactions.
                'dust_emit'   The photons were last emitted dust and did
                              not undergo any subsequent interactions
                'source_scat' The photons were last emitted from a
                              source and were subsequently scattered
                'dust_scat'   The photons were last emitted from dust and
                              were subsequently scattered

        inclination : int, optional
            The number of the viewing angle to plot (zero-based). Use 'all'
            to return all viewing angles.

        vmin, vmax : float, optional
            The range in fluxes to show (in ergs/s if distance is not
            specified, or in ergs/cm^2/s if distance is specified).

        cmap : matplotlib colormap
            The color to plot the SED(s) in.

        labels : bool, optional
            Whether or not to show the axis labels

        source_id, dust_id : int, optional
            If the output file was made with track_origin='detailed', a
            specific source and dust component can be specified. If
            these are not specified, then the total component requested
            for all sources or dust types is returned.

        Returns
        -------

        axes : matplotlib.pyplot.Axes instance (if `axes` was set)
            The updated matplotlib Axes are returned if they were passed in
            as input.

        figure : matplotlib.pyplot.Figure instance (if `axes` and `filename` were not set)
            The matplotlib Figure created
        '''

        if axes is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = axes

        wav, nufnu = self.get_image(stokes=stokes, group=group, technique=technique, distance=distance, component=component, inclination=inclination, uncertainties=False)

        iw = np.argmin(np.abs(wav - wavelength))

        ax.imshow(nufnu[iw, :, :], vmin=vmin, vmax=vmax, cmap=cmap)

        if filename:
            fig.savefig(filename)
            plt.close(fig)
            return
        elif axes is None:
            return fig
        else:
            return ax

    def get_available_components(self, iteration=-1):
        '''
        Find out what physical components are available in the output file

        Parameters
        ----------

        iteration : integer, optional
            The iteration to retrieve the grid for. The default is to return the components for the last iteration
        '''

        # Open file
        f = h5py.File('%s.rtout' % self.name, 'r')

        # If iteration is last one, find iteration number
        if iteration == -1:
            iteration = find_last_iteration(f)

        # Return components
        components = f['Iteration %05i' % iteration].keys()
        if 'specific_energy' in components:
            components.append('temperature')
        return components

    def get_physical_grid(self, name, iteration=-1, dust_id='all'):
        '''
        Retrieve one of the physical grids for the model

        Parameters
        ----------

        name : str
            The component to retrieve. This should be one of:
                'specific_energy'  The specific energy absorbed in each cell
                'temperature'      The dust temperature in each cell (only
                                   available for cells with LTE dust)
                'density'          The density in each cell (after possible
                                   dust sublimation)
                'density_diff'     The difference in the final density
                                   compared to the initial density
                'n_photons'        The number of unique photons that went
                                   through each cell

        iteration : integer, optional
            The iteration to retrieve the grid for. The default is to return the grid for the last iteration.

        Returns
        -------

        array : numpy.array instance for regular grids
            The physical grid in cgs.

        Notes
        -----

        At the moment, this method only works on regular grids, not AMR or Oct-tree grids
        '''

        # Check name
        available_components = self.get_available_components()
        if name not in available_components:
            raise Exception("name should be one of %s" % string.join(available_components, '/'))

        # Open file
        f_out = h5py.File('%s.rtout' % self.name, 'r')

        # If iteration is last one, find iteration number
        if iteration == -1:
            iteration = find_last_iteration(f_out)

        # Extract specific energy grid
        if name == 'temperature':
            array = np.array(f_out['Iteration %05i' % iteration]['specific_energy'])
            f_in = h5py.File('%s.rtin' % self.name, 'r')
            g_dust = f_in['Dust']
            for i in range(array.shape[0]):
                dust = g_dust['dust_%03i' % (i + 1)].file
                d = SphericalDust(dust)
                array[i, :, :, :] = d.mean_opacities._specific_energy2temperature(array[i, :, :, :])
        else:
            array = np.array(f_out['Iteration %05i' % iteration][name])

        # If required, extract grid for a specific dust type
        if dust_id == 'all':
            return [array[i, :, :, :] for i in range(array.shape[0])]
        else:
            return array[dust_id, :, :, :]
