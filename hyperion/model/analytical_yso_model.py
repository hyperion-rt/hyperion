from __future__ import print_function, division

import collections
from copy import deepcopy

import numpy as np
from astropy import log as logger

from ..densities import FlaredDisk, AlphaDisk, PowerLawEnvelope, UlrichEnvelope, AmbientMedium
from ..util.interpolate import interp1d_fast_loglog
from ..util.constants import pi, sigma, c, G
from ..sources import SphericalSource, SpotSource
from ..util.functions import FreezableClass, virtual_file
from ..grid import SphericalPolarGrid, CylindricalPolarGrid

from . import Model


def _min_none(*args):
    "Minimum of several arguments, ignoring None values"
    return min(x for x in args if x is not None)


def _max_none(*args):
    "Maximum of several arguments, ignoring None values"
    return max(x for x in args if x is not None)


class Star(FreezableClass):

    def __init__(self):
        self.sources = {}
        self.sources['star'] = SphericalSource(name='star')
        self.mass = None
        self.radius = None
        self.limb = False
        self._freeze()

    @classmethod
    def read(self, filename):
        raise Exception("Can only call ``read`` for Model, not AnalyticalYSOModel")

    def add_spot(self, *args, **kwargs):
        self.sources['star'].add_spot(SpotSource(*args, **kwargs))

    def __setattr__(self, attribute, value):
        if self.isfinal():
            raise Exception("Attribute %s can no longer be changed" % attribute)
        if attribute in ['luminosity', 'temperature', 'spectrum']:
            self.sources['star'].__setattr__(attribute, value)
            return
        elif attribute in ['radius', 'limb']:
            for source in self.sources:
                self.sources[source].__setattr__(attribute, value)
        FreezableClass.__setattr__(self, attribute, value)

    def __getattr__(self, attribute):
        if attribute in ['luminosity', 'temperature', 'spectrum', 'radius', 'limb']:
            return getattr(self.sources['star'], attribute)
        else:
            return FreezableClass.__getattr__(self, attribute)

    def total_luminosity(self):
        "Return the total luminosity of the star, including accretion"
        ltot = 0.
        for source in self.sources:
            if self.sources[source].luminosity is not None:
                ltot += self.sources[source].luminosity
        return ltot

    def effective_temperature(self):
        "Return the effective temperature of the star, including accretion"
        return (self.total_luminosity() / (4. * pi * self.radius ** 2. * sigma)) ** 0.25

    def total_spectrum(self, bnu_range=None):
        '''
        Return the total spectrum of the star, including accretion

        Parameters
        ----------
        bnu_range : tuple
            Range of frequencies to cover for sources that have Planck spectra
        '''

        # Retrieve all the spectra for the sources of emission
        nu_all, fnu_all = [], []
        for source in self.sources:
            if self.sources[source].temperature is not None:
                if bnu_range is None:
                    raise ValueError("bnu_range is needed for sources with Planck spectra")
                nu, fnu = self.sources[source].get_spectrum(nu_range=bnu_range)
            else:
                nu, fnu = self.sources[source].get_spectrum()
            nu_all.append(nu)
            fnu_all.append(fnu)

        # Find common minimum and maximum for all spectra
        nu_min = np.min([nu.min() for nu in nu_all])
        nu_max = np.max([nu.max() for nu in nu_all])

        # Find common frequencies
        nu_common = np.unique(np.sort(np.hstack(nu_all)))
        nu_common = nu_common[(nu_common >= nu_min) & (nu_common <= nu_max)]

        # Compute total spectrum
        fnu_total = np.zeros(nu_common.shape)
        for i in range(len(self.sources)):
            fnu_total += interp1d_fast_loglog(nu_all[i], fnu_all[i], nu_common, bounds_error=False, fill_value=0.)

        return nu_common, fnu_total


class AnalyticalYSOModel(Model):

    def __init__(self, name=None):
        "Initialize an analytical YSO model"

        self.star = Star()
        self.disks = []
        self.envelopes = []
        self.ambients = []

        Model.__init__(self, name=name)

    def add_density_grid(self, *args, **kwargs):
        raise NotImplementedError("add_density_grid cannot be used for AnalyticalYSOModel")

    def use_quantities(self, *args, **kwargs):
        raise NotImplementedError("use_quantities cannot be used for AnalyticalYSOModel")

    def use_geometry(self, *args, **kwargs):
        raise NotImplementedError("use_geometry cannot be used for AnalyticalYSOModel")

    # DENSITY COMPONENTS

    def add_ambient_medium(self, subtract=[]):
        '''
        Add an ambient density medium to the model

        Parameters
        ----------
        subtract : list
            Components to subtract from the ambient density medium (see
            notes below).

        Returns
        -------
        ambient : :class:`~hyperion.densities.AmbientMedium`
            An :class:`~hyperion.densities.AmbientMedium` instance.

        Examples
        --------

        To add an ambient medium to the model, you can do::

            >>> ambient = m.add_ambient_medium()

        then set the ambient medium properties using e.g.::

            >>> from hyperion.util.constants import au, pc
            >>> ambient.rho = 1.e-20  # cgs
            >>> ambient.rmin = 0.1 * au  # cm
            >>> ambient.rmax = pc  # cm

        See the :class:`~hyperion.densities.AmbientMedium` documentation
        to see which parameters can be set.

        Notes
        -----

        By default, the ambient medium simply adds a constant density ``rho`` of
        dust to the whole model between the inner and outer radius. However, it
        is possible to pass components that should be subtracted from the
        constant density using the ``subtract=`` argument. In the following
        example::

            >>> e = m.add_power_law_envelope()
            >>> m.add_ambient_medium(subtract=[e])

        the ambient medium does not simply add a constant density ``rho`` of
        dust everywhere, but it adds dust such that the density never falls
        below ``rho`` between ``rmin`` and ``rmax`` - that is, it subtracts the
        density of component ``e`` from the ``rho``, with a minimum of zero. In
        regions where the density of component of ``e`` is larger than ``rho``,
        no dust is added.
        '''
        ambient = AmbientMedium()
        ambient.star = self.star
        ambient.subtract = subtract
        self.ambients.append(ambient)
        return ambient

    def add_flared_disk(self):
        '''
        Add a flared disk to the model

        Returns
        -------
        disk : :class:`~hyperion.densities.FlaredDisk`
            A :class:`~hyperion.densities.FlaredDisk` instance.

        Examples
        --------

        To add a flared disk to the model, you can do::

            >>> disk = m.add_flared_disk()

        then set the disk properties using e.g.::

            >>> disk.mass = 1.e30  # g
            >>> disk.rmin = 1e10  # cm
            >>> disk.rmax = 1e14  # cm

        See the :class:`~hyperion.densities.FlaredDisk` documentation
        to see which parameters can be set.
        '''
        disk = FlaredDisk()
        disk.star = self.star
        self.disks.append(disk)
        return disk

    def add_alpha_disk(self):
        '''
        Add an alpha disk to the geometry

        This is similar to a flared disk, but with accretion luminosity. See
        :class:`~hyperion.densities.AlphaDisk` for more details.

        Returns
        -------
        disk : :class:`~hyperion.densities.AlphaDisk`
            A :class:`~hyperion.densities.AlphaDisk` instance.

        Examples
        --------

        To add an alpha disk to the model, you can do::

            >>> disk = m.add_alpha_disk()

        then set the disk properties using e.g.::

            >>> disk.mass = 1.e30  # g
            >>> disk.rmin = 1e10  # cm
            >>> disk.rmax = 1e14  # cm

        See the :class:`~hyperion.densities.AlphaDisk` documentation
        to see which parameters can be set.
        '''
        disk = AlphaDisk()
        disk.star = self.star
        self.disks.append(disk)
        return disk

    def add_settled_disks(self, reference_disk, reference_size, eta=0.,
                          sizes=[], dust_files=[]):
        '''
        Automatically create disks with varying degrees of settling

        .. warning:: this function is still experimental, and will be documented once stable
        '''

        exists = False

        for disk in self.disks:
            if disk is reference_disk:
                logger.warn("Reference disk already exists, not re-adding")
                exists = True

        if not exists:
            logger.warn("Reference disk does not exist, adding")
            self.disks.append(reference_disk)

        for i, size in enumerate(sizes):
            disk = deepcopy(reference_disk)
            disk.star = self.star
            disk.h_0 *= (size / reference_size) ** (-eta)
            disk.dust = dust_files[i]
            self.disks.append(disk)

    def add_ulrich_envelope(self):
        '''
        Add an infalling rotationally flatted envelope to the model

        Returns
        -------
        env : :class:`~hyperion.densities.UlrichEnvelope`
            An :class:`~hyperion.densities.UlrichEnvelope` instance.

        Examples
        --------

        To add an infalling envelope to the model, you can do::

            >>> env = m.add_ulrich_envelope()

        then set the envelope properties using e.g.::

            >>> from hyperion.util.constants import msun, yr, au
            >>> env.mdot = 1.e-6 * msun / yr  # g/s
            >>> env.rmin = 0.1 * au  # cm
            >>> env.rmax = 10000. * au  # cm

        See the :class:`~hyperion.densities.UlrichEnvelope` documentation
        to see which parameters can be set.
        '''
        envelope = UlrichEnvelope()
        envelope.star = self.star
        self.envelopes.append(envelope)
        return envelope

    def add_power_law_envelope(self):
        '''
        Add a spherically symmetric power-law envelope to the model

        Returns
        -------
        env : :class:`~hyperion.densities.PowerLawEnvelope`
            A :class:`~hyperion.densities.PowerLawEnvelope` instance.

        Examples
        --------

        To add a power-law envelope to the model, you can do::

            >>> env = m.add_power_law_envelope()

        then set the envelope properties using e.g.::

            >>> from hyperion.util.constants import msun, au
            >>> env.mass = 0.1 * msun  # g/s
            >>> env.rmin = 0.1 * au  # cm
            >>> env.rmax = 10000. * au  # cm

        See the :class:`~hyperion.densities.PowerLawEnvelope` documentation
        to see which parameters can be set.
        '''
        envelope = PowerLawEnvelope()
        envelope.star = self.star
        self.envelopes.append(envelope)
        return envelope

    def _check_all_set(self):
        for disk in self.disks:
            disk._check_all_set()
        for envelope in self.envelopes:
            envelope._check_all_set()
        for ambient in self.ambients:
            ambient._check_all_set()

    # MIDPLANE OPTICAL DEPTH

    def print_midplane_tau(self, wavelength):
        for i, disk in enumerate(self.disks):
            if disk.mass > 0.:
                tau = (disk.midplane_cumulative_density(np.array([disk.rmax]))
                       * disk.dust.interp_chi_wav(wavelength))
                print("Disk %i: %.5e" % (i + 1, tau))

    def get_midplane_tau(self, r):

        self._check_all_set()

        # Find the total combined optical depth through the midplane
        tau_midplane = np.zeros(r.shape)

        # IMPLEMENT: PAHs

        for i, disk in enumerate(self.disks):
            if disk.mass > 0.:
                if disk.dust is None:
                    raise Exception("Disk %i dust not set" % i)
                nu_min = disk.dust.optical_properties.nu[0]
                nu_max = disk.dust.optical_properties.nu[-1]
                nu, fnu = self.star.total_spectrum(bnu_range=[nu_min, nu_max])
                if np.any(fnu > 0.):
                    tau_midplane += (disk.midplane_cumulative_density(r)
                                     * disk.dust.optical_properties.chi_planck_spectrum(nu, fnu))

        for i, envelope in enumerate(self.envelopes):
            if envelope.exists():
                if envelope.dust is None:
                    raise Exception("envelope %i dust not set" % i)
                nu_min = envelope.dust.optical_properties.nu[0]
                nu_max = envelope.dust.optical_properties.nu[-1]
                nu, fnu = self.star.total_spectrum(bnu_range=[nu_min, nu_max])
                if np.any(fnu > 0.):
                    tau_midplane += (envelope.midplane_cumulative_density(r)
                                     * envelope.dust.optical_properties.chi_planck_spectrum(nu, fnu))

        return tau_midplane

    def plot_midplane_tau(self, filename):

        import matplotlib.pyplot as plt

        tau_midplane = self.get_midplane_tau(self.grid.r_wall[1:])

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.loglog(self.grid.r[1:] / self.grid.r[1] - 1.,
                  tau_midplane[1:] - tau_midplane[:-1],
                  drawstyle='steps-mid')
        fig.savefig(filename)

    # COORDINATE GRID

    def radial_range(self):

        if (len(self.disks) == 0 and
            len(self.envelopes) == 0 and
            len(self.ambients) == 0):
            rmin = self.star.radius
        else:
            rmin_values = ([disk.rmin for disk in self.disks] +
                           [envelope.rmin for envelope in self.envelopes] +
                           [ambient.rmin for ambient in self.ambients])
            rmin = _min_none(*rmin_values)

        rmax_values = [self.star.radius]
        rmax_values += ([disk.rmax for disk in self.disks] +
                        [envelope.rmax for envelope in self.envelopes] +
                        [ambient.rmax for ambient in self.ambients])
        rmax = _max_none(*rmax_values)

        return rmin, rmax

    def set_cylindrical_polar_grid_auto(self, n_w, n_z, n_phi,
                                        wmax=None, zmax=None, min_spacing=1.e-8):
        '''
        Set the grid to be cylindrical polar with automated resolution.

        Parameters
        ----------
        n_w, n_z, n_phi : int
            Number of cells to use in the radial, vertical, and azimuthal
            directions.
        wmax : float, optional
            The maximum radius to extend out to. If not specified, this is
            set to the maximum  cylindrical radius of the dust geometry in the
            mid-plane.
        zmax : float, optional
            The maximum height above and below the midplane to extend to. If
            not specified, this is set to the maximum cylindrical radius of
            the dust geometry.
        min_spacing : float, optional
            The minimum spacing (in relative terms) for the inner cell walls.
            The spacing from rmin to the next cell wall cannot be smaller than
            rmin * (1 + min_spacing).
        '''
        self.grid = {'grid_type': 'cylindrical',
                     'n1': n_w, 'n2': n_z, 'n3': n_phi,
                     'rmax': wmax, 'zmax': zmax, 'min_spacing':min_spacing}

    def set_spherical_polar_grid_auto(self, n_r, n_theta, n_phi,
                                      rmax=None, min_spacing=1.e-8):
        '''
        Set the grid to be spherical polar with automated resolution.

        Parameters
        ----------
        n_r, n_theta, n_phi : int
            Number of cells to use in the radial, theta, and azimuthal
            directions.
        rmax : float, optional
            The maximum radius to extend out to. If not specified, this is
            set to the maximum spherical radius of the dust geometry in the
            mid-plane. Note that if you are including a disk with a
            cylindrical outer edge, this should be set to a value larger
            than the disk radius, otherwise the disk will be truncated with
            a spherical edge.
        min_spacing : float, optional
            The minimum spacing (in relative terms) for the inner cell walls.
            The spacing from rmin to the next cell wall cannot be smaller than
            rmin * (1 + min_spacing).
        '''
        self.grid = {'grid_type': 'spherical',
                     'n1': n_r, 'n2': n_theta, 'n3': n_phi,
                     'rmax': rmax, 'min_spacing':min_spacing}

    def _set_polar_grid_auto(self, n1=None, n2=None, n3=None, grid_type=None,
                             zmax=None, rmax=None, min_spacing=1.e-8):

        if self.star.radius is None:
            raise Exception("The central source radius need to be defined "
                            "before the grid can be set up")

        if grid_type is 'spherical':
            n_r, n_theta, n_phi = n1, n2, n3
        elif grid_type is 'cylindrical':
            n_r, n_z, n_phi = n1, n2, n3
        else:
            raise Exception("Unknown grid type: %s" % grid_type)

        # Find minimum and maximum radius
        if len(self.disks) == 0 and len(self.envelopes) == 0:
            rmin = self.star.radius
        else:
            rmin_values = ([disk.rmin for disk in self.disks] +
                           [envelope.rmin for envelope in self.envelopes] +
                           [ambient.rmin for ambient in self.ambients])
            rmin = _min_none(*rmin_values)

        if not rmax:
            rmax_values = [2. * self.star.radius]
            rmax_values += ([disk.rmax for disk in self.disks] +
                            [envelope.rmax for envelope in self.envelopes] +
                            [ambient.rmax for ambient in self.ambients])
            rmax = _max_none(*rmax_values)

        if rmax < rmin:
            logger.warn("Grid rmax < rmin, model will consist only of central star")
            rmin = self.star.radius
            rmax = 2. * self.star.radius

        if np.isnan(rmin):
            raise Exception("R_min is NaN")

        if np.isnan(rmax):
            raise Exception("R_max is NaN")

        if rmin == 0:
            raise ValueError("R_min is 0, so cannot set up the grid cell "
                             "walls automatically. Use set_%s_polar_grid()"
                             " instead to specify the cell wall positions"
                             "." % grid_type)

        # RADIAL WALLS

        # Set first wall to be at half the stellar radius to avoid any
        # numerical problems. Set second wall to be the disk inner radius.

        r_wall = [self.star.radius / 2., rmin]

        # Define a radial grid to compute the midplane column density on
        r = np.logspace(-20., np.log10((rmax - rmin) / rmin), 100000) * rmin + rmin

        # We need the first point to be exactly at the inner radius
        r[0] = rmin

        # Get cumulative midplane optical depth
        tau_midplane = self.get_midplane_tau(r)

        # Find where the second wall after rmin would be if grid cells were
        # spaced equally
        r_next_real = rmin * ((rmax / rmin) ** (1. / n_r) - 1.)

        # Find where the second wall would be if we put it at tau=0.1
        if tau_midplane[-1] <= 0.1:
            r_next_tau = rmax - rmin
        else:
            r_next_tau = np.interp(0.1, tau_midplane, r) - rmin

        # Pick the smallest
        rnext = min(r_next_real, r_next_tau)

        # Make sure rnext isn't too small
        if rmin * (1. + min_spacing) > rnext + rmin:
            logger.warn("Spacing of inner radial cells is too small, resetting to {0}".format(min_spacing))
            rnext = rmin * min_spacing

        # Define wall positions
        r_wall = np.hstack([0., np.logspace(np.log10(rnext / rmin), np.log10((rmax - rmin) / rmin), n_r - 1)]) * rmin + rmin
        r_wall = np.hstack([0., r_wall])

        if grid_type is 'spherical':

            # THETA WALLS
            t_wall = np.linspace(0, pi, n_theta + 1)
            t_wall = t_wall + np.sin(2 * t_wall) / 6.

        elif grid_type is 'cylindrical':

            # Z WALLS
            zmin = np.inf
            for disk in self.disks:
                zmin = min(zmin, disk.scale_height_at(rmin))

            if not zmax:
                zmax = rmax
            if n_z % 2 == 0:
                n_zn = n_z / 2
                z_wall1 = np.linspace(zmin * 0.1, zmin * 0.9, 10)
                z_wall2 = np.logspace(np.log10(zmin),
                                      np.log10(zmax),
                                      n_zn - 10)
                z_wall = np.hstack([z_wall1, z_wall2])
                z_wall = np.hstack([-z_wall[::-1], z_wall])
            else:
                n_zn = (n_z - 1) / 2
                z_wall1 = np.linspace(zmin * 0.1, zmin * 0.9, 10)
                z_wall2 = np.logspace(np.log10(zmin),
                                      np.log10(zmax),
                                      n_zn - 10)
                z_wall = np.hstack([z_wall1, z_wall2])
                z_wall = np.hstack([-z_wall[::-1], 0., z_wall])

        # PHI WALLS
        p_wall = np.linspace(0., 2. * pi, n_phi + 1)

        if grid_type is 'spherical':
            return SphericalPolarGrid(r_wall, t_wall, p_wall)
        else:
            return CylindricalPolarGrid(r_wall, z_wall, p_wall)

    # ACCRETION

    def setup_magnetospheric_accretion(self, mdot, rtrunc, fspot,
                                       xwav_min=0.001, xwav_max=0.01):
        '''
        Set up the model for magnetospheric accretion

        Parameters
        ----------
        mdot : float
            The accretion rate onto the star in cgs
        rtrunc : float
            The magnetospheric truncation radius of the disk in cgs
        fspot : float
            The spot coverage fraction. Photons will be emitted uniformly from
            the star, the coverage fraction ``fspot`` will determine the
            spectrum of the hot spot emission (smaller covering fractions will
            lead to a hotter spectrum).

        Notes
        -----
        This method only takes into account the hot spot and X-ray emission
        from the stellar surface. To simulate the viscous accretion luminosity
        in the disk, add an :class:`~hyperion.densities.AlphaDisk` to the
        model using :meth:`~hyperion.model.AnalyticalYSOModel.add_alpha_disk`
        and set the accretion rate or luminosity accordingly.
        '''

        # For convenience
        lstar = self.star.sources['star'].luminosity

        if self.star.mass is None:
            raise Exception("Stellar mass is not set")

        # Find the luminosity dissipated in the shock
        lshock = G * self.star.mass * mdot * (1 / self.star.radius - 1 / rtrunc)

        # Hot spot parameters
        fluxratio = 0.5 * lshock / lstar / fspot
        teff = (lstar / (4. * pi * self.star.radius ** 2 * sigma)) ** 0.25  # Kelvin
        tshock = teff * (1 + fluxratio) ** 0.25  # Kelvin

        # Set the hot spot source
        self.star.sources['uv'] = SphericalSource(name='uv',
                                                  radius=self.star.radius)
        self.star.sources['uv'].luminosity = lshock / 2. + lstar * fspot
        self.star.sources['uv'].temperature = tshock

        # X-rays from 0.1 to 10nm
        wav = np.logspace(np.log10(xwav_min), np.log10(xwav_max), 100)[::-1]
        nu = c * 1.e4 / wav
        fnu = np.repeat(1., nu.shape)

        # Set the X-ray source
        self.star.sources['xray'] = SphericalSource(name='xray',
                                                    radius=self.star.radius)
        self.star.sources['xray'].luminosity = lshock / 2.
        self.star.sources['xray'].spectrum = (nu, fnu)

        # Reduce the total luminosity from the original source
        self.star.sources['star'].luminosity *= 1 - fspot

    # OUTPUT

    def to_model(self, merge_if_possible=True):
        '''
        Returns a Model instance of the current model

        The AnalyticalYSOModel class is dynamic in the sense that one can
        change the parameters relating to the density structure at any time.
        This method computes the Model instance corresponding to the current
        density structure and computes the optimal grid.

        Parameters
        ----------
        merge_if_possible : bool
            Whether to merge density arrays that have the same dust type
        '''

        if self.grid is None:
            raise Exception("The coordinate grid needs to be defined")

        # Initialize new Model instance
        m = Model()

        # Set up grid
        if isinstance(self.grid, collections.Mapping):
            m.grid = self._set_polar_grid_auto(**self.grid)
        else:
            m.grid = deepcopy(self.grid)

        # Copy over settings
        m.name = self.name
        m.conf = deepcopy(self.conf)
        m.sources = deepcopy(self.sources)
        m.binned_output = deepcopy(self.binned_output)
        m.peeled_output = deepcopy(self.peeled_output)

        m._minimum_temperature = deepcopy(self._minimum_temperature)
        m._minimum_specific_energy = deepcopy(self._minimum_specific_energy)

        m._monochromatic = self._monochromatic
        m._frequencies = self._frequencies

        # Easiest way to copy over all settings
        g = virtual_file()
        self.write_run_conf(g)
        m.read_run_conf(g)
        g.close()

        for i, disk in enumerate(self.disks):

            if disk.rmin >= disk.rmax:
                logger.warn("Disk rmin >= rmax, ignoring density contribution")
            elif disk.mass == 0.:
                logger.warn("Disk mass is zero, ignoring density contribution")
            else:

                if not disk.dust:
                    raise Exception("Disk %i dust not set" % (i + 1))
                m.add_density_grid(disk.density(m.grid), disk.dust,
                                   merge_if_possible=merge_if_possible)

        for i, envelope in enumerate(self.envelopes):

            if envelope.rmin >= envelope.rmax:
                logger.warn("Envelope rmin >= rmax, "
                            "ignoring density contribution")
            elif isinstance(envelope, UlrichEnvelope) and envelope.rho_0 == 0.:
                logger.warn("Ulrich envelope has zero density everywhere, "
                            "ignoring density contribution")
            elif isinstance(envelope, PowerLawEnvelope) and envelope.mass == 0.:
                logger.warn("Power-law envelope has zero density everywhere, "
                            "ignoring density contribution")
            else:

                if not envelope.dust:
                    raise Exception("Envelope dust not set")
                m.add_density_grid(envelope.density(m.grid), envelope.dust,
                                   merge_if_possible=merge_if_possible)

                if envelope.cavity is not None:
                    if envelope.cavity.theta_0 == 0.:
                        logger.warn("Cavity opening angle is zero, "
                                    "ignoring density contribution")
                    elif envelope.cavity.rho_0 == 0.:
                        logger.warn("Cavity density is zero, "
                                    "ignoring density contribution")
                    else:
                        if not envelope.cavity.dust:
                            raise Exception("Cavity dust not set")
                        m.add_density_grid(envelope.cavity.density(m.grid),
                                           envelope.cavity.dust,
                                           merge_if_possible=merge_if_possible)

        # AMBIENT MEDIUM

        for i, ambient in enumerate(self.ambients):

            if ambient.density == 0.:
                logger.warn("Ambient medium has zero density, "
                            "ignoring density contribution")
            else:

                if not ambient.dust:
                    raise Exception("Ambient medium dust not set")

                m.add_density_grid(ambient.density(m.grid),
                                   ambient.dust,
                                   merge_if_possible=merge_if_possible)

        # SOURCES

        # Star

        if self.star.sources['star'].luminosity > 0:
            if self.star.sources['star'] not in self.sources:
                m.add_source(self.star.sources['star'])

        # Accretion

        if 'uv' in self.star.sources and self.star.sources['uv'].luminosity > 0.:
            if self.star.sources['uv'] not in self.sources:
                m.add_source(self.star.sources['uv'])

        if 'xray' in self.star.sources and self.star.sources['xray'].luminosity > 0.:
            if self.star.sources['xray'] not in self.sources:
                m.add_source(self.star.sources['xray'])

        for i, disk in enumerate(self.disks):

            if isinstance(disk, AlphaDisk):
                if disk.rmin >= disk.rmax:
                    logger.warn("Disk rmin >= rmax, "
                                "ignoring accretion luminosity")
                elif disk.mass == 0.:
                    logger.warn("Disk mass is zero, "
                                "ignoring accretion luminosity")
                elif disk.lvisc == 0.:
                    logger.warn("Disk viscous luminosity is zero, "
                                "ignoring accretion luminosity")
                else:
                    m.add_map_source(luminosity=disk.lvisc, map=disk.accretion_luminosity(m.grid), name='accdisk%i' % i)

        return m

    def write(self, filename=None, compression=True, copy=True,
              absolute_paths=False, wall_dtype=float,
              physics_dtype=float, overwrite=True, merge_if_possible=True):
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
        merge_if_possible : bool
            Whether to merge density arrays that have the same dust type
        '''

        self.filename = filename

        m = self.to_model(merge_if_possible=merge_if_possible)

        m.write(filename=filename, compression=compression,
                copy=copy, absolute_paths=absolute_paths,
                wall_dtype=wall_dtype, physics_dtype=physics_dtype,
                overwrite=overwrite)
