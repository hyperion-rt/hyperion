from __future__ import print_function, division

import numpy as np

from ..util.functions import FreezableClass, bool2str, str2bool, is_numpy_array


class OutputConf(FreezableClass):

    def __init__(self):
        '''
        Initialize default output configuration
        '''
        self.output_density = 'none'
        self.output_density_diff = 'none'
        self.output_specific_energy = 'last'
        self.output_n_photons = 'none'
        self._freeze()

    @classmethod
    def read(cls, group):
        self = cls()
        self.output_density = group.attrs['output_density'].decode('utf-8')
        self.output_density_diff = group.attrs['output_density_diff'].decode('utf-8')
        self.output_specific_energy = group.attrs['output_specific_energy'].decode('utf-8')
        self.output_n_photons = group.attrs['output_n_photons'].decode('utf-8')
        return self

    def write(self, group):
        group.attrs['output_density'] = np.string_(self.output_density.encode('utf-8'))
        group.attrs['output_density_diff'] = np.string_(self.output_density_diff.encode('utf-8'))
        group.attrs['output_specific_energy'] = np.string_(self.output_specific_energy.encode('utf-8'))
        group.attrs['output_n_photons'] = np.string_(self.output_n_photons.encode('utf-8'))


class RunConf(object):

    def __init__(self):
        '''
        Initialize default run configuration
        '''
        self.set_n_initial_iterations(5)
        self.n_photons = {}
        self.set_propagation_check_frequency(0.001)
        self.set_seed(-124902)
        self.set_raytracing(False)
        self.set_max_interactions(1000000)
        self.set_max_reabsorptions(1000000)
        self.set_pda(False)
        self.set_mrw(False)
        self.set_convergence(False)
        self.set_kill_on_absorb(False)
        self.set_forced_first_scattering(True)
        self.set_output_bytes(8)
        self.set_sample_sources_evenly(False)
        self.set_enforce_energy_range(True)
        self.set_copy_input(True)
        self._monochromatic = False
        super(RunConf, self).__init__()

    def set_propagation_check_frequency(self, frequency):
        '''
        Set how often to check that the photon is in the right cell

        During photon propagation, it is possible that floating point issues
        cause a photon to end up in the wrong cell. By default, the code will
        randomly double check the position and cell of a photon for every 1 in
        1000 cell wall crossings, but this can be adjusted with this method.
        Note that values higher than around 0.001 will cause the code to slow
        down.

        Parameters
        ----------
        frequency : float
            How often the photon position and cell should be double-checked (1
            is always, 0 is never).
        '''
        if not np.isscalar(frequency) or isinstance(frequency, basestring):
            raise TypeError("frequency should be a scalar value")
        if frequency < 0. or frequency > 1.:
            raise ValueError("frequency should be between 0 and 1")
        self._frequency = frequency

    def _read_propagation_check_frequency(self, group):
        self._frequency = group.attrs['propagation_check_frequency']

    def _write_propagation_check_frequency(self, group):
        group.attrs['propagation_check_frequency'] = self._frequency

    def set_seed(self, seed):
        '''
        Set the seed for the random number generation

        Parameters
        ----------
        seed : int
            The seed with which to initialize the random number generation.
            This should be negative.
        '''
        if type(seed) != int or seed >= 0:
            raise ValueError("seed should be a negative integer")
        self._seed = seed

    def _read_seed(self, group):
        self._seed = group.attrs['seed']

    def _write_seed(self, group):
        group.attrs['seed'] = self._seed

    def set_n_initial_iterations(self, n_iter):
        '''
        Set the number of initial iterations for computing the specific
        energy in each cell.

        Parameters
        ----------
        n_iter : int
            The number of initial iterations
        '''
        self.n_iter = n_iter

    def _read_n_initial_iterations(self, group):
        self.n_iter = group.attrs['n_initial_iter']

    def _write_n_initial_iterations(self, group):
        group.attrs['n_initial_iter'] = self.n_iter

    def set_n_photons(self, initial=None, imaging=None,
                      imaging_sources=None, imaging_dust=None,
                      raytracing_sources=None, raytracing_dust=None,
                      stats=0):
        '''
        Set the number of photons for the different iterations

        Note that any values not specified will be set to the default
        values.

        Parameters
        ----------
        initial : float, optional
            Number of photons for the initial specific energy iterations
        imaging : float, optional
            Number of photons for the main SED/image iteration. This argument
            is used in the case of non-monochromatic radiation transfer.
        imaging_sources : float, optional
            Number of photons emitted from sources during the main SED/image
            iteration in the case of monochromatic radiation transfer.
        imaging_dust : float, optional
            Number of photons emitted from dust during the main SED/image
            iteration in the case of monochromatic radiation transfer.
        raytracing_sources : float, optional
            Number of photons emitted from sources during the raytracing
            SED/image iteration, if applicable.
        raytracing_dust : float, optional
            Number of photons emitted from dust during the raytracing
            SED/image iteration, if applicable.
        stats : float, optional
            How often to print out statistics. Also used to determine the
            photon chunk size for MPI.
        '''

        if self.n_iter == 0:
            if initial is not None:
                raise Exception("[n_photons] initial should not be set since no initial interations are being computed")
            else:
                if 'initial' in self.n_photons:
                    del self.n_photons['initial']
        else:
            if initial is None:
                raise Exception("[n_photons] initial should be set since the initial iterations are being computed")
            else:
                self.n_photons['initial'] = initial

        if self.raytracing:
            if raytracing_sources is None:
                raise Exception("[n_photons] raytracing_sources needs to be set in raytracing mode")
            else:
                self.n_photons['raytracing_sources'] = raytracing_sources
            if raytracing_dust is None:
                raise Exception("[n_photons] raytracing_dust needs to be set in raytracing mode")
            else:
                self.n_photons['raytracing_dust'] = raytracing_dust
        else:
            if raytracing_sources is not None:
                raise Exception("[n_photons] raytracing_sources should not be set as raytracing is not being used")
            if raytracing_dust is not None:
                raise Exception("[n_photons] raytracing_dust should not be set as raytracing is not being used")

        if self._monochromatic:
            if imaging_sources is None:
                raise Exception("[n_photons] imaging_sources needs to be set in monochromatic mode")
            else:
                self.n_photons['last_sources'] = imaging_sources
            if imaging_dust is None:
                raise Exception("[n_photons] imaging_dust needs to be set in monochromatic mode")
            else:
                self.n_photons['last_dust'] = imaging_dust
            if imaging is not None:
                raise Exception("[n_photons] imaging should not be set in monochromatic mode")
        else:
            if imaging_sources is not None:
                raise Exception("[n_photons] imaging_sources should not be set as the monochromatic option is not being used")
            if imaging_dust is not None:
                raise Exception("[n_photons] imaging_dust should not be set as the monochromatic option is not being used")
            if imaging is None:
                raise Exception("[n_photons] imaging should bet set")
            else:
                self.n_photons['last'] = imaging

        self.n_photons['stats'] = stats

    def _read_n_photons(self, group):

        if self.n_iter != 0:
            self.n_photons['initial'] = group.attrs['n_initial_photons']

        if self._monochromatic:
            self.n_photons['last_sources'] = group.attrs['n_last_photons_sources']
            self.n_photons['last_dust'] = group.attrs['n_last_photons_dust']
        else:
            self.n_photons['last'] = group.attrs['n_last_photons']

        if self.raytracing:
            self.n_photons['raytracing_sources'] = group.attrs['n_ray_photons_sources']
            self.n_photons['raytracing_dust'] = group.attrs['n_ray_photons_dust']

        self.n_photons['stats'] = group.attrs['n_stats']

    def _write_n_photons(self, group):

        if self.n_photons == {}:
            raise Exception("Photon numbers not set")

        if self.n_iter == 0:
            if 'initial' in self.n_photons and self.n_photons['initial'] is not None:
                raise Exception("[n_photons] initial should not be set since no initial interations are being computed")
        else:
            if 'initial' in self.n_photons and self.n_photons['initial'] is not None:
                group.attrs['n_initial_photons'] = self.n_photons['initial']
            else:
                raise Exception("[n_photons] initial should be set since the initial iterations are being computed")

        if self._monochromatic:
            if 'last_sources' in self.n_photons:
                group.attrs['n_last_photons_sources'] = self.n_photons['last_sources']
            else:
                raise Exception("[n_photons] imaging_sources needs to be set in monochromatic mode")
            if 'last_dust' in self.n_photons:
                group.attrs['n_last_photons_dust'] = self.n_photons['last_dust']
            else:
                raise Exception("[n_photons] imaging_dust needs to be set in monochromatic mode")
            if 'last' in self.n_photons:
                raise Exception("[n_photons] imaging should not be set in monochromatic mode")
        else:
            if 'last_sources' in self.n_photons:
                raise Exception("[n_photons] imaging_sources should not be set as the monochromatic option is not being used")
            if 'last_dust' in self.n_photons:
                raise Exception("[n_photons] imaging_dust should not be set as the monochromatic option is not being used")
            if 'last' in self.n_photons:
                group.attrs['n_last_photons'] = self.n_photons['last']
            else:
                raise Exception("[n_photons] imaging should bet set")

        if self.raytracing:
            if 'raytracing_sources' in self.n_photons:
                group.attrs['n_ray_photons_sources'] = self.n_photons['raytracing_sources']
            else:
                raise Exception("[n_photons] raytracing_sources needs to be set in raytracing mode")
            if 'raytracing_dust' in self.n_photons:
                group.attrs['n_ray_photons_dust'] = self.n_photons['raytracing_dust']
            else:
                raise Exception("[n_photons] raytracing_dust needs to be set in raytracing mode")
        else:
            if 'raytracing_sources' in self.n_photons:
                raise Exception("[n_photons] raytracing_sources should not be set as raytracing is not being used")
            if 'raytracing_dust' in self.n_photons:
                raise Exception("[n_photons] raytracing_dust should not be set as raytracing is not being used")

        group.attrs['n_stats'] = self.n_photons['stats']

    def set_raytracing(self, raytracing):
        '''
        Set whether to use raytracing for the non-scattered flux

        If enabled, only scattered photons are peeled off in the iteration
        following the initial iterations, and an additional final
        iteration is carrried out, with raytracing of the remaining flux
        (sources and thermal and non-thermal dust emission).

        Parameters
        ----------
        raytracing : bool
            Whether or not to use raytracing in the final iteration
        '''
        self.raytracing = raytracing

    def _read_raytracing(self, group):
        self.raytracing = str2bool(group.attrs['raytracing'])

    def _write_raytracing(self, group):
        group.attrs['raytracing'] = bool2str(self.raytracing)

    def set_max_interactions(self, inter_max):
        '''
        Set the maximum number of interactions a photon can have.

        Parameters
        ----------
        inter_max : int
            Maximum number of interactions for a single photon. This can be
            used to prevent photons from getting stuck in very optically
            thick regions, especially if the modified random walk is not
            used.
        '''
        self.n_inter_max = inter_max

    def _read_max_interactions(self, group):
        self.n_inter_max = group.attrs['n_inter_max']

    def _write_max_interactions(self, group):
        group.attrs['n_inter_max'] = self.n_inter_max

    def set_max_reabsorptions(self, reabs_max):
        '''
        Set the maximum number of successive reabsorptions by a source that a
        photon can have.

        Parameters
        ----------
        reabs_max : int
            Maximum number of reabsorptions for a single photon.
        '''
        self.n_reabs_max = reabs_max

    def _read_max_reabsorptions(self, group):
        self.n_reabs_max = group.attrs['n_reabs_max']

    def _write_max_reabsorptions(self, group):
        group.attrs['n_reabs_max'] = self.n_reabs_max

    def set_pda(self, pda):
        '''
        Set whether to use the Partial Diffusion Approximation (PDA)

        If enabled, the PDA is used to compute the specific energy in cells
        which have seen few or no photons by formally solving the diffusion
        equations, using the cells with valid specific energies as boundary
        conditions.

        Parameters
        ----------
        pda : bool
            Whether or not to use the PDA

        References
        ----------
        Min et al. 2009, Astronomy and Astrophysics, 497, 155
        '''
        self.pda = pda

    def _read_pda(self, group):
        self.pda = str2bool(group.attrs['pda'])

    def _write_pda(self, group):
        group.attrs['pda'] = bool2str(self.pda)

    def set_mrw(self, mrw, gamma=1.0, inter_max=1000):
        '''
        Set whether to use the Modified Random Walk (MRW) approximation

        If enabled, the MRW speeds up the propagation of photons in very
        optically thick regions by locally setting up a spherical diffusion
        region.

        Parameters
        ----------
        mrw : bool
            Whether or not to use the MRW
        gamma : float, optional
            The parameter describing the starting criterion for the MRW.
            The MRW is carried out if the distance to the closest cell is
            larger than `gamma` times the Rosseland mean free path.
        inter_max : int, optional
            Maximum number of interactions during a single random walk.
            This can be used to prevent photons from getting stuck in the
            corners of cells in very optically thick regions, where the MRW
            stars to become inefficient itself.

        References
        ----------
        Min et al. 2009, Astronomy and Astrophysics, 497, 155
        '''
        self.mrw = mrw
        self.mrw_gamma = gamma
        self.n_inter_mrw_max = inter_max

    def _read_mrw(self, group):
        self.mrw = str2bool(group.attrs['mrw'])
        if self.mrw:
            self.mrw_gamma = group.attrs['mrw_gamma']
            self.n_inter_mrw_max = group.attrs['n_inter_mrw_max']

    def _write_mrw(self, group):
        group.attrs['mrw'] = bool2str(self.mrw)
        if self.mrw:
            group.attrs['mrw_gamma'] = self.mrw_gamma
            group.attrs['n_inter_mrw_max'] = self.n_inter_mrw_max

    def set_convergence(self, convergence, percentile=100., absolute=0., relative=0.):
        '''
        Set whether to check for convergence over the initial iterations

        If enabled, the code will check whether the specific energy absorbed
        in each cell has converged. First, the ratio between the previous
        and current specific energy absorbed in each cell is computed in each
        cell, and the value at the specified percentile (`percentile`) is
        found. Then, convergence has been achieved if this value is less than
        an absolute threshold (`absolute`), and if it changed by less than
        a relative threshold ratio (`relative`).

        Parameters
        ----------
        convergence : bool
            Whether or not to check for convergence.
        percentile : float, optional
            The percentile at which to check for convergence.
        absolute : float, optional
            The abolute threshold below which the percentile value of the
            ratio has to be for convergence.
        relative : float, optional
            The relative threshold below which the ratio in the percentile
            value has to be for convergence.
        '''
        self.check_convergence = True
        self.convergence_percentile = percentile
        self.convergence_absolute = absolute
        self.convergence_relative = relative

    def _read_convergence(self, group):
        self.check_convergence = str2bool(group.attrs['check_convergence'])
        if self.check_convergence:
            self.convergence_percentile = group.attrs['convergence_percentile']
            self.convergence_absolute = group.attrs['convergence_absolute']
            self.convergence_relative = group.attrs['convergence_relative']

    def _write_convergence(self, group):
        group.attrs['check_convergence'] = bool2str(self.check_convergence)
        if self.check_convergence:
            group.attrs['convergence_percentile'] = self.convergence_percentile
            group.attrs['convergence_absolute'] = self.convergence_absolute
            group.attrs['convergence_relative'] = self.convergence_relative

    def set_kill_on_absorb(self, kill_on_absorb):
        '''
        Set whether to kill absorbed photons

        Parameters
        ----------
        kill_on_absorb : bool
            Whether to kill absorbed photons
        '''
        self.kill_on_absorb = kill_on_absorb

    def _read_kill_on_absorb(self, group):
        self.kill_on_absorb = str2bool(group.attrs['kill_on_absorb'])

    def _write_kill_on_absorb(self, group):
        group.attrs['kill_on_absorb'] = bool2str(self.kill_on_absorb)

    def set_forced_first_scattering(self, forced_first_scattering):
        '''
        Set whether to ensure that photons scatter at least once before
        escaping the grid.

        Parameters
        ----------
        forced_first_scattering : bool
            Whether to force at least one scattering before escaping the
            grid

        References
        ----------
        Wood & Reynolds, 1999, The Astrophysical Journal, 525, 799
        '''
        self.forced_first_scattering = forced_first_scattering

    def _read_forced_first_scattering(self, group):
        self.forced_first_scattering = str2bool(group.attrs['forced_first_scattering'])

    def _write_forced_first_scattering(self, group):
        group.attrs['forced_first_scattering'] = bool2str(self.forced_first_scattering)

    def set_enforce_energy_range(self, enforce):
        '''
        Set how to deal with cells that have specific energy rates that are
        below or above that provided in the mean opacities and emissivities.

        Parameters
        ----------
        enforce : bool
            Whether or not to reset specific energies that are above or below
            the range of values used to specify the mean opacities and
            emissivities to the maximum or minimum value of the range. Setting
            this to True modifies the energy in the simulation, but ensures
            that the emissivities are consistent with the energy in the cells.
            Setting this to False means that the total energy in the grid will
            be correct, but that the emissivities may be inconsistent with the
            energy in the cells (if an energy is out of range, the code will
            pick the closest available one). In both cases, warnings will be
            displayed to notify the user whether this is happening.
        '''

        self.enforce_energy_range = enforce

    def _read_enforce_energy_range(self, group):
        self.enforce_energy_range = str2bool(group.attrs['enforce_energy_range'])

    def _write_enforce_energy_range(self, group):
        group.attrs['enforce_energy_range'] = bool2str(self.enforce_energy_range)

    def set_copy_input(self, copy):
        '''
        Set whether to copy the input data into the output file.

        Parameters
        ----------
        copy : bool
            Whether to copy the input data into the output file (True) or
            whether to link to it (False)
        '''
        self.copy_input = copy

    def _read_copy_input(self, group):
        self.copy_input = str2bool(group.attrs['copy_input'])

    def _write_copy_input(self, group):
        group.attrs['copy_input'] = bool2str(self.copy_input)

    def set_output_bytes(self, io_bytes):
        '''
        Set whether to output physical quantity arrays in 32-bit or 64-bit

        Parameters
        ----------
        io_bytes : int
            The number of bytes for the output. This should be either 4
            (for 32-bit) or 8 (for 64-bit).
        '''
        self.physics_io_bytes = io_bytes

    def _read_output_bytes(self, group):
        self.physics_io_bytes = group.attrs['physics_io_bytes']

    def _write_output_bytes(self, group):
        group.attrs['physics_io_bytes'] = self.physics_io_bytes

    def set_sample_sources_evenly(self, sample_sources_evenly):
        '''
        If set to 'True', sample evenly from all sources and apply
        probability weight based on relative luminosities. Otherwise,
        sample equal energy photons from sources with probability given by
        relative luminosities.

        Parameters
        ----------
        sample_evenly : bool
            Whether to sample different sources evenly
        '''
        self.sample_sources_evenly = sample_sources_evenly

    def _read_sample_sources_evenly(self, group):
        self.sample_sources_evenly = str2bool(group.attrs['sample_sources_evenly'])

    def _write_sample_sources_evenly(self, group):
        group.attrs['sample_sources_evenly'] = bool2str(self.sample_sources_evenly)

    def read_run_conf(self, group):  # not a class method because inherited
        '''
        Read the configuation in from an HDF5 group

        Parameters
        ----------
        group : h5py.highlevel.File or h5py.highlevel.Group
            The HDF5 group to read the configuration from
        '''
        self._read_propagation_check_frequency(group)
        self._read_seed(group)
        self._read_n_initial_iterations(group)
        self._read_raytracing(group)
        self._read_n_photons(group)
        self._read_max_interactions(group)
        self._read_max_reabsorptions(group)
        self._read_pda(group)
        self._read_mrw(group)
        self._read_convergence(group)
        self._read_kill_on_absorb(group)
        self._read_forced_first_scattering(group)
        self._read_output_bytes(group)
        self._read_sample_sources_evenly(group)
        self._read_enforce_energy_range(group)
        self._read_copy_input(group)

    def write_run_conf(self, group):
        '''
        Writes out the configuation to an HDF5 group

        Parameters
        ----------
        group : h5py.highlevel.File or h5py.highlevel.Group
            The HDF5 group to write the configuration to
        '''
        self._write_propagation_check_frequency(group)
        self._write_seed(group)
        self._write_n_initial_iterations(group)
        self._write_n_photons(group)
        self._write_raytracing(group)
        self._write_max_interactions(group)
        self._write_max_reabsorptions(group)
        self._write_pda(group)
        self._write_mrw(group)
        self._write_convergence(group)
        self._write_kill_on_absorb(group)
        self._write_forced_first_scattering(group)
        self._write_output_bytes(group)
        self._write_sample_sources_evenly(group)
        self._write_enforce_energy_range(group)
        self._write_copy_input(group)


class ImageConf(FreezableClass):

    def __init__(self, sed=True, image=True):
        '''
        Initialize default image configuration
        '''
        self.sed = sed
        self.image = image
        if self.sed:
            self.set_aperture_range(1, np.inf, np.inf)  # reasonable default
        # Since there are no reasonable defaults for image size and limits,
        # as well as wavelength range, we don't set any. But we still have
        # to set the variables since new attributes cannot be created later.
        if self.image:
            self.n_x = self.n_y = None
            self.xmin = self.xmax = self.ymin = self.ymax = None
        self.n_wav = self.wav_min = self.wav_max = None
        self.set_output_bytes(8)
        self.set_track_origin('no')
        self.set_uncertainties(False)
        self._monochromatic = False
        self._freeze()

    def set_output_bytes(self, io_bytes):
        '''
        Set whether to output images/SEDs in 32-bit or 64-bit.

        Parameters
        ----------
        io_bytes : int
            The number of bytes for the output. This should be either 4
            (for 32-bit) or 8 (for 64-bit).
        '''
        self.io_bytes = io_bytes

    def _read_output_bytes(self, group):
        self.io_bytes = group.attrs['io_bytes']

    def _write_output_bytes(self, group):
        group.attrs['io_bytes'] = self.io_bytes

    def set_image_size(self, n_x, n_y):
        '''
        Set the size of the output images

        Parameters
        ----------
        n_x, n_y : int
            The number of pixels in the x and y directions
        '''
        if type(n_x) != int:
            raise Exception("n_x should be an integer")
        if type(n_y) != int:
            raise Exception("n_y should be an integer")
        self.n_x = n_x
        self.n_y = n_y

    def _read_image_size(self, group):
        self.n_x = group.attrs['n_x']
        self.n_y = group.attrs['n_y']

    def _write_image_size(self, group):
        if self.n_x is None or self.n_y is None:
            raise Exception("Image size has not been set")
        group.attrs['n_x'] = self.n_x
        group.attrs['n_y'] = self.n_y

    def set_image_limits(self, xmin, xmax, ymin, ymax):
        '''
        Set the extent of the output images

        Parameters
        ----------
        xmin, xmax, ymin, ymax : float
            The extent of the images, which are either in cm (if using
            standard binned images or peeloff images) or in degrees (if
            using peeling off to a point inside the model).
        '''
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def _read_image_limits(self, group):
        self.xmin = group.attrs['x_min']
        self.xmax = group.attrs['x_max']
        self.ymin = group.attrs['y_min']
        self.ymax = group.attrs['y_max']

    def _write_image_limits(self, group):
        if self.xmin is None or self.xmax is None or self.ymin is None or self.ymax is None:
            raise Exception("Image limits have not been set")
        group.attrs['x_min'] = self.xmin
        group.attrs['x_max'] = self.xmax
        group.attrs['y_min'] = self.ymin
        group.attrs['y_max'] = self.ymax

    def set_aperture_range(self, n_ap, ap_min, ap_max):
        '''
        Set the range of apertures to use for SEDs

        Parameters
        ----------
        n_ap : int
            The number of apertures to compute SEDs in
        ap_min, ap_max : float
            The smallest and largest aperture to use, in cm
        '''
        if type(n_ap) is not int:
            raise Exception("n_ap should be an integer")
        self.n_ap = n_ap
        self.ap_min = ap_min
        self.ap_max = ap_max

    def _read_aperture_range(self, group):
        self.n_ap = group.attrs['n_ap']
        self.ap_min = group.attrs['ap_min']
        self.ap_max = group.attrs['ap_max']

    def _write_aperture_range(self, group):
        group.attrs['n_ap'] = self.n_ap
        group.attrs['ap_min'] = self.ap_min
        group.attrs['ap_max'] = self.ap_max

    def set_wavelength_range(self, n_wav, wav_min, wav_max):
        '''
        Set the range of wavelengths to use for SEDs

        Parameters
        ----------
        n_wav : int
            The number of wavelengths to compute SEDs for
        wav_min, wav_max : float
            The smallest and largest wavelength to use, in microns
        '''
        if type(n_wav) is not int:
            raise Exception("n_wav should be an integer")
        self.n_wav = n_wav
        self.wav_min = wav_min
        self.wav_max = wav_max

    def _read_wavelength_range(self, group):

        self.n_wav = group.attrs['n_wav']

        if 'inu_min' in group.attrs:
            self.wav_min = group.attrs['inu_min']
            self.wav_max = group.attrs['inu_max']
        else:
            self.wav_min = group.attrs['wav_min']
            self.wav_max = group.attrs['wav_max']

    def _write_wavelength_range(self, group):

        if self.n_wav is None:
            raise Exception("Wavelength range has not been set")
        group.attrs['n_wav'] = self.n_wav

        if self._monochromatic:
            if not type(self.wav_min) is int or not type(self.wav_max) is int:
                raise Exception("In monochromatic mode, wavelength range should be given as integers")
            group.attrs['inu_min'] = self.wav_min
            group.attrs['inu_max'] = self.wav_max
        else:
            group.attrs['wav_min'] = self.wav_min
            group.attrs['wav_max'] = self.wav_max

    def set_track_origin(self, track_origin):
        '''
        Set whether to track the origin of the photons. The options are:

        'no' - does not split up the images/SEDs by origin

        'basic' - this splits up the images/SEDs into:

          * The photons last emitted from a source and did not undergo
            any subsequent interactions.
          * The photons last emitted from dust and did not undergo any
            subsequent interactions
          * The photons last emitted from a source and subsequently
            scattered
          * The photons last emitted from dust and subsequently
            scattered

        'detailed' - as above, but in each category, the origin is further
                     refined into each individual source and dust type.

        Parameters
        ----------
        track_origin : bool
            Whether to track the origin of the photons as described above.
        '''

        if track_origin is True:
            track_origin = 'basic'
        elif track_origin is False:
            track_origin = 'no'
        elif track_origin not in ['no', 'basic', 'detailed']:
            raise Exception("track_origin should be one of no/basic/detailed")

        self.track_origin = track_origin

    def _read_track_origin(self, group):
        self.track_origin = group.attrs['track_origin'].decode('ascii')

    def _write_track_origin(self, group):
        group.attrs['track_origin'] = np.string_(self.track_origin.encode('utf-8'))

    def set_uncertainties(self, uncertainties):
        '''
        Set whether to compute uncertainties on the images/SEDs

        Parameters
        ----------
        uncertainties : bool
            Whether to compute uncertainties on the images/SEDs.
        '''
        self.uncertainties = uncertainties

    def _read_uncertainties(self, group):
        self.uncertainties = str2bool(group.attrs['uncertainties'])

    def _write_uncertainties(self, group):
        group.attrs['uncertainties'] = bool2str(self.uncertainties)

    @classmethod
    def read(cls, group):
        self = cls()
        self._read_viewing_info(group)
        self._read_main_info(group)
        return self

    def write(self, group):
        self._write_viewing_info(group)
        self._write_main_info(group)

    def _read_viewing_info(self, group):
        pass

    def _write_viewing_info(self, group):
        pass

    def _read_main_info(self, group):
        self.sed = str2bool(group.attrs['compute_sed'])
        self.image = str2bool(group.attrs['compute_image'])
        if self.sed:
            self._read_aperture_range(group)
        if self.image:
            self._read_image_size(group)
            self._read_image_limits(group)
        self._read_wavelength_range(group)
        self._read_output_bytes(group)
        self._read_track_origin(group)
        self._read_uncertainties(group)

    def _write_main_info(self, group):
        group.attrs['compute_sed'] = bool2str(self.sed)
        group.attrs['compute_image'] = bool2str(self.image)
        if self.sed:
            self._write_aperture_range(group)
        if self.image:
            self._write_image_size(group)
            self._write_image_limits(group)
        self._write_wavelength_range(group)
        self._write_output_bytes(group)
        self._write_track_origin(group)
        self._write_uncertainties(group)


class BinnedImageConf(ImageConf):

    def __init__(self, n_theta=None, n_phi=None, **kwargs):
        self.n_theta = n_theta
        self.n_phi = n_phi
        ImageConf.__init__(self, **kwargs)

    def set_viewing_bins(self, n_theta, n_phi):
        '''
        Set the number of viewing angles to use

        Parameters
        ----------
        n_theta, n_phi
            The number of viewing angles to use in the theta and phi
            directions respectively.
        '''
        self.n_theta = n_theta
        self.n_phi = n_phi

    def _read_viewing_bins(self, group):
        self.n_theta = group.attrs['n_theta']
        self.n_phi = group.attrs['n_phi']

    def _write_viewing_bins(self, group):
        group.attrs['n_theta'] = self.n_theta
        group.attrs['n_phi'] = self.n_phi

    def _read_viewing_info(self, group):
        self._read_viewing_bins(group)

    def _write_viewing_info(self, group):
        self._write_viewing_bins(group)


class PeeledImageConf(ImageConf):

    def __init__(self, **kwargs):
        self.set_viewing_angles([], [])
        self.set_inside_observer(None)
        self.set_ignore_optical_depth(False)
        self.set_peeloff_origin(None)
        self.set_depth(-np.inf, np.inf)
        ImageConf.__init__(self, **kwargs)

    def set_viewing_angles(self, theta, phi):
        '''
        Set the viewing angles to use

        Parameters
        ----------
        theta, phi : iterable of floats
            The viewing angles to compute SEDs for.

        Examples
        -------

        Set viewing angles using lists of well-defined angles:

        >>> image.set_viewing_angles([30.,55.,87.],[22.,44.,34.])

        Set viewing angles using generated numpy arrays:

        >>> image.set_viewing_angles(np.linspace(0.,90.,10), np.repeat(30.,10))

        Set a single viewing direction:

        >>> image.set_viewing_angles([77.],[25.])
        '''

        if isinstance(theta, (list, tuple)):
            theta = np.array(theta)
        if isinstance(phi, (list, tuple)):
            phi = np.array(phi)

        if not is_numpy_array(theta) or theta.ndim != 1:
            raise ValueError("theta should be a 1-D sequence")
        if not is_numpy_array(phi) or phi.ndim != 1:
            raise ValueError("phi should be a 1-D sequence")

        if len(theta) != len(phi):
            raise ValueError("Length of theta and phi arrays do not match")

        self.viewing_angles = list(zip(theta, phi))
        self.n_view = len(self.viewing_angles)

    def _read_viewing_angles(self, group):
        angles = group['angles']
        self.viewing_angles = list(zip(angles['theta'], angles['phi']))

    def _write_viewing_angles(self, group):
        group.attrs['n_view'] = len(self.viewing_angles)
        group.create_dataset('angles', data=np.array(self.viewing_angles, dtype=[('theta', float), ('phi', float)]))

    def set_inside_observer(self, position):
        '''
        Set the observer to be inside the model

        Parameters
        ----------
        position : tuple of 3 floats
           The spatial coordinates of the observer, in cm
        '''
        if isinstance(position, (list, tuple)):
            position = np.array(position)
        if position is not None:
            if not is_numpy_array(position) or position.ndim != 1 or len(position) != 3:
                raise ValueError("position should be a 1-D sequence with 3 elements")
        self.inside_observer = position

    def _read_inside_observer(self, group):
        self.inside_observer = (group.attrs['observer_x'],
                                group.attrs['observer_y'],
                                group.attrs['observer_z'])

    def _write_inside_observer(self, group):
        group.attrs['observer_x'] = self.inside_observer[0]
        group.attrs['observer_y'] = self.inside_observer[1]
        group.attrs['observer_z'] = self.inside_observer[2]

    def set_ignore_optical_depth(self, ignore_optical_depth):
        '''
        Ingore optical depth when creating images.

        This is useful in cases where one wants to understand how much the
        optical depth is affecting a set of images.

        Parameters
        ----------
        ignore_optical_depth : bool
           Whether to ignore optical depth effects (default is False)
        '''
        self.ignore_optical_depth = ignore_optical_depth

    def _read_ignore_optical_depth(self, group):
        self.ignore_optical_depth = str2bool(group.attrs['ignore_optical_depth'])

    def _write_ignore_optical_depth(self, group):
        group.attrs['ignore_optical_depth'] = bool2str(self.ignore_optical_depth)

    def set_peeloff_origin(self, position):
        '''
        Set the origin for the peeloff.

        Parameters
        ----------
        position : tuple of 3 floats
           The coordinates of the origin of the peeling-off, in cm
        '''
        if isinstance(position, (list, tuple)):
            position = np.array(position)
        if position is not None:
            if not is_numpy_array(position) or position.ndim != 1 or len(position) != 3:
                raise ValueError("position should be a 1-D sequence with 3 elements")
        self.peeloff_origin = position

    def _read_peeloff_origin(self, group):
        self.peeloff_origin = (group.attrs['peeloff_x'],
                               group.attrs['peeloff_y'],
                               group.attrs['peeloff_z'])

    def _write_peeloff_origin(self, group):
        group.attrs['peeloff_x'] = self.peeloff_origin[0]
        group.attrs['peeloff_y'] = self.peeloff_origin[1]
        group.attrs['peeloff_z'] = self.peeloff_origin[2]

    def set_depth(self, d_min, d_max):
        '''
        Set the minimum and maximum distance between which photons should be
        peeled off.

        By default, d_min and d_max are set to -inf and +inf respectively.
        This option can be useful to compute for example models in a spherical
        polar grid, but including only the photons in a slab.

        Parameters
        ----------
        d_min, d_max : float
           The minimum and maximum distance between which photons should be
           peeled-off. Distance increases away from the observer, and d_min
           and d_max are the distances closest and furthest from the observer
           respectively. The origin is the position of the observer if inside
           the model, otherwise it is the origin of the grid.
        '''
        self.d_min = d_min
        self.d_max = d_max

    def _read_depth(self, group):
        self.d_min = group.attrs['d_min']
        self.d_max = group.attrs['d_max']

    def _write_depth(self, group):
        group.attrs['d_min'] = self.d_min
        group.attrs['d_max'] = self.d_max

    def _read_viewing_info(self, group):

        if str2bool(group.attrs['inside_observer']):
            self._read_inside_observer(group)
        else:
            self._read_peeloff_origin(group)

        self._read_ignore_optical_depth(group)

        self._read_viewing_angles(group)

        self._read_depth(group)

    def _write_viewing_info(self, group):

        if self.peeloff_origin is not None and self.inside_observer is not None:
            raise Exception("Cannot specify inside observer and peeloff origin at the same time")

        if self.inside_observer is not None:
            group.attrs['inside_observer'] = bool2str(True)
            self._write_inside_observer(group)
            if self.viewing_angles == []:
                self.set_viewing_angles([90.], [0.])
            if self.image and self.xmin < self.xmax:
                raise ValueError("longitudes should increase towards the left for inside observers")
            if self.d_min < 0.:
                if self.d_min != -np.inf:
                    raise ValueError("Lower limit of depth should be positive for inside observer")
                self.d_min = 0.
            if self.d_max < 0.:
                raise ValueError("Upper limit of depth should be positive for inside observer")

        elif len(self.viewing_angles) > 0:
            group.attrs['inside_observer'] = bool2str(False)
            if self.peeloff_origin is None:
                self.set_peeloff_origin((0., 0., 0.))
            self._write_peeloff_origin(group)
        else:
            raise Exception("Need to specify either observer position, or viewing angles")

        self._write_ignore_optical_depth(group)

        self._write_viewing_angles(group)

        self._write_depth(group)
