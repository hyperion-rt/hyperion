import numpy as np
from hyperion.util.functions import FreezableClass


def bool2str(value):
    if value:
        return "yes"
    else:
        return "no"


class OutputConf(FreezableClass):

    def __init__(self):
        '''
        Initialize default output configuration
        '''
        self.output_temperature = 'all'
        self.output_density = 'all'
        self.output_density_diff = 'all'
        self.output_specific_energy_abs = 'all'
        self.output_n_photons = 'all'
        self._freeze()

    def write(self, group):
        group.attrs['output_temperature'] = self.output_temperature
        group.attrs['output_density'] = self.output_density
        group.attrs['output_density_diff'] = self.output_density_diff
        group.attrs['output_specific_energy_abs'] = self.output_specific_energy_abs
        group.attrs['output_n_photons'] = self.output_n_photons


class RunConf(FreezableClass):

    def __init__(self):
        '''
        Initialize default run configuration
        '''
        self.set_n_temperature_iterations(5)
        self.set_n_photons()
        self.set_raytracing(False)
        self.set_max_interactions(1000000)
        self.set_max_reabsorptions(1000000)
        self.set_pda(False)
        self.set_mrw(False)
        self.set_convergence(False)
        self.set_kill_on_absorb(False)
        self.set_forced_first_scattering(True)
        self.set_dust_sublimation('no')
        self.set_minimum_temperature(0.1)
        self.set_output_bytes(8)
        self._freeze()

    def set_n_temperature_iterations(self, n_iter):
        '''
        Set the number of temperature iterations

        Parameters
        ----------
        n_iter : int
            The number of temperature iterations
        '''
        self.n_iter = n_iter

    def _write_n_temperature_iterations(self, group):
        group.attrs['n_lucy_iter'] = self.n_iter

    def set_n_photons(self, temperature=None, imaging=None,
                      imaging_mono_sources=None, imaging_mono_dust=None,
                      raytracing_sources=None, raytracing_dust=None, stats=10000):
        '''
        Set the number of photons for the different iterations

        Note that any values not specified will be set to the default
        values.

        Parameters
        ----------
        temperature : float, optional
            Number of photons for the temperature iterations
        imaging : float, optional
            Number of photons for the main SED/image iteration. This argument
            is used in the case of non-monochromatic radiation transfer.
        imaging_mono_sources : float, optional
            Number of photons emitted from sources during the main SED/image
            iteration in the case of monochromatic radiation transfer.
        imaging_mono_dust : float, optional
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

        if imaging is not None and imaging_mono_sources is not None:
            raise Exception("imaging and imaging_mono_sources cannot both be set")

        if imaging is not None and imaging_mono_dust is not None:
            raise Exception("imaging and imaging_mono_dust cannot both be set")

        self.n_lucy_photons = temperature
        self.n_last_photons = imaging
        self.n_last_photons_sources = imaging_mono_sources
        self.n_last_photons_dust = imaging_mono_dust
        self.n_ray_photons_sources = raytracing_sources
        self.n_ray_photons_dust = raytracing_dust
        self.n_stats = stats

    def _write_n_photons(self, group):

        if self.n_lucy_photons is not None:
            group.attrs['n_lucy_photons'] = self.n_lucy_photons
        if self.n_last_photons is not None:
            group.attrs['n_last_photons'] = self.n_last_photons
        if self.n_last_photons_sources is not None:
            group.attrs['n_last_photons_sources'] = self.n_last_photons_sources
        if self.n_last_photons_dust is not None:
            group.attrs['n_last_photons_dust'] = self.n_last_photons_dust
        if self.n_ray_photons_sources is not None:
            group.attrs['n_ray_photons_sources'] = self.n_ray_photons_sources
        if self.n_ray_photons_dust is not None:
            group.attrs['n_ray_photons_dust'] = self.n_ray_photons_dust
        if self.n_stats is not None:
            group.attrs['n_stats'] = self.n_stats

    def set_raytracing(self, raytracing):
        '''
        Set whether to use raytracing for the non-scattered flux

        If enabled, only scattered photons are peeled off in the iteration
        following the temperature calculation, and an additional final
        iteration is carrried out, with raytracing of the remaining flux
        (sources and thermal and non-thermal dust emission).

        Parameters
        ----------
        raytracing : bool
            Whether or not to use raytracing in the final iteration
        '''
        self.raytracing = raytracing

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

    def _write_max_reabsorptions(self, group):
        group.attrs['n_reabs_max'] = self.n_reabs_max

    def set_pda(self, pda):
        '''
        Set whether to use the Partial Diffusion Approximation (PDA)

        If enabled, the PDA is used to compute the temperature in cells
        which have seen few or no photons by formally solving the diffusion
        equations, using the cells with valid temperatures as boundary
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

    def _write_mrw(self, group):
        group.attrs['mrw'] = bool2str(self.mrw)
        if(self.mrw):
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

    def _write_convergence(self, group):
        group.attrs['check_convergence'] = bool2str(self.check_convergence)
        if(self.check_convergence):
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
        self.kill_on_absorb = bool2str(kill_on_absorb)

    def _write_kill_on_absorb(self, group):
        group.attrs['kill_on_absorb'] = self.kill_on_absorb

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
        self.forced_first_scattering = bool2str(forced_first_scattering)

    def _write_forced_first_scattering(self, group):
        group.attrs['forced_first_scattering'] = self.forced_first_scattering

    def set_dust_sublimation(self, mode, temperature=1600):
        '''
        Set the dust sublimation parameters.

        Parameters
        ----------
        mode : str
            The dust sublimation mode, which can be:
                * 'no'   - no sublimation
                * 'fast' - remove all dust in cells exceeding the
                           sublimation temperature
                * 'slow' - reduce the dust in cells exceeding the
                           sublimation temperature
        temperature : float, optional
            The dust sublimation temperature, in K
        '''
        self.dust_sublimation_mode = mode
        self.dust_sublimation_temperature = temperature

    def _write_dust_sublimation(self, group):
        group.attrs['dust_sublimation_mode'] = self.dust_sublimation_mode
        if self.dust_sublimation_mode in ['slow', 'fast', 'cap']:
            group.attrs['dust_sublimation_temperature'] = self.dust_sublimation_temperature

    def set_minimum_temperature(self, temperature):
        '''
        Set the minimum dust temperature

        Dust which has a temperature that falls below this value will be
        reset to the minimum at the end of each iteration.

        Parameters
        ----------
        temperature : float
            The minimum temperature in K
        '''
        self.minimum_temperature = temperature

    def _write_minimum_temperature(self, group):
        group.attrs['minimum_temperature'] = self.minimum_temperature

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

    def _write_output_bytes(self, group):
        group.attrs['physics_io_bytes'] = self.physics_io_bytes

    def write(self, group):
        '''
        Writes out the configuation to an HDF5 group

        Parameters
        ----------
        group : h5py.highlevel.File or h5py.highlevel.Group
            The HDF5 group to write the configuration to
        '''
        self._write_n_temperature_iterations(group)
        self._write_n_photons(group)
        self._write_raytracing(group)
        self._write_max_interactions(group)
        self._write_max_reabsorptions(group)
        self._write_pda(group)
        self._write_mrw(group)
        self._write_convergence(group)
        self._write_kill_on_absorb(group)
        self._write_forced_first_scattering(group)
        self._write_dust_sublimation(group)
        self._write_minimum_temperature(group)
        self._write_output_bytes(group)


class ImageConf(FreezableClass):

    def __init__(self, sed=True, image=True):
        '''
        Initialize default image configuration
        '''
        self.sed = sed
        self.image = image
        if self.sed:
            self.set_aperture_range(1, np.inf, np.inf)
        if self.image:
            self.set_image_size(1, 1)
            self.set_image_limits(-np.inf, np.inf, -np.inf, np.inf)
        self.set_wavelength_range(250, 0.01, 5000.)
        self.set_output_bytes(8)
        self.set_track_origin(False)
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
        self.n_x = n_x
        self.n_y = n_y

    def _write_image_size(self, group):
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

    def _write_image_limits(self, group):
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

    def _write_wavelength_range(self, group):
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
        Set whether to track the origin of the photons. This splits them up
        into:

          * The photons were last emitted from a source and did not undergo
            any subsequent interactions.
          * The photons were last emitted dust and did not undergo any
            subsequent interactions
          * The photons were last emitted from a source and were
            subsequently scattered
          * The photons were last emitted from dust and were subsequently
            scattered

        Parameters
        ----------
        track_origin : bool
            Whether to track the origin of the photons as described above.
        '''
        self.track_origin = track_origin

    def _write_track_origin(self, group):
        group.attrs['track_origin'] = bool2str(self.track_origin)

    def set_uncertainties(self, uncertainties):
        '''
        Set whether to compute uncertainties on the images/SEDs

        Parameters
        ----------
        uncertainties : bool
            Whether to compute uncertainties on the images/SEDs.
        '''
        self.uncertainties = uncertainties

    def _write_uncertainties(self, group):
        group.attrs['uncertainties'] = bool2str(self.uncertainties)

    def write(self, group):
        self._write_viewing_info(group)
        self._write_main_info(group)

    def _write_viewing_info(self, group):
        pass

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

    def _write_viewing_bins(self, group):
        group.attrs['n_theta'] = self.n_theta
        group.attrs['n_phi'] = self.n_phi

    def _write_viewing_info(self, group):
        self._write_viewing_bins(group)


class PeeledImageConf(ImageConf):

    def __init__(self, viewing_angles=None, inside_observer=None, **kwargs):
        if viewing_angles is not None:
            self.n_view = len(self.viewing_angles)
        else:
            self.n_view = 0
        self.viewing_angles = viewing_angles
        self.inside_observer = inside_observer
        self.set_peeloff_origin((0., 0., 0.))
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
        if len(theta) != len(phi):
            raise Exception("Length of theta and phi arrays do not match")
        self.viewing_angles = zip(theta, phi)
        self.n_view = len(self.viewing_angles)

    def _write_viewing_angles(self, group):
        group.attrs['n_view'] = len(self.viewing_angles)
        group.create_dataset('Angles', data=np.array(self.viewing_angles, dtype=[('theta', float), ('phi', float)]))

    def set_inside_observer(self, position):
        '''
        Set the observer to be inside the model

        Parameters
        ----------
        position : tuple of 3 floats
           The coordinates of the observer, in cm
        '''
        self.inside_observer = position

    def _write_inside_observer(self, group):
        group.attrs['observer_x'] = self.inside_observer[0]
        group.attrs['observer_y'] = self.inside_observer[1]
        group.attrs['observer_z'] = self.inside_observer[2]

    def set_peeloff_origin(self, position):
        '''
        Set the origin for the peeloff.

        Parameters
        ----------
        position : tuple of 3 floats
           The coordinates of the origin of the peeling-off, in cm
        '''
        self.peeloff_origin = position

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

    def _write_depth(self, group):
        group.attrs['d_min'] = self.d_min
        group.attrs['d_max'] = self.d_max

    def _write_viewing_info(self, group):

        if self.viewing_angles and self.inside_observer:
            raise Exception("Cannot specify inside observer and viewing angles at the same time")

        if self.inside_observer is not None:
            group.attrs['inside_observer'] = 'yes'
            self._write_inside_observer(group)
        elif self.viewing_angles is not None:
            group.attrs['inside_observer'] = 'no'
            self._write_viewing_angles(group)
            self._write_peeloff_origin(group)
        else:
            raise Exception("Need to specify either observer position, or viewing angles")

        self._write_depth(group)
