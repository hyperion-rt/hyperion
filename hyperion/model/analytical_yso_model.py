from copy import deepcopy

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as mpl

from hyperion.model import Model
from hyperion.densities import FlaredDisk, AlphaDiskWhitney, PowerLawEnvelope, UlrichEnvelope, AmbientMedium
from hyperion.util.interpolate import interp1d_fast_loglog
from hyperion.util.constants import pi, sigma, c, G
from hyperion.sources import SphericalSource, SpotSource
from hyperion.util.functions import FreezableClass
from hyperion.util.convenience import OptThinRadius
from hyperion.util.logger import logger


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

    def total_spectrum(self):
        "Return the total spectrum of the star, including accretion"

        # Retrieve all the spectra for the sources of emission
        nu_all, fnu_all = [], []
        for source in self.sources:
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

    def _finalize(self):
        for source in self.sources:
            self.sources[source]._finalize()
        FreezableClass._finalize(self)


class AnalyticalYSOModel(Model):

    def __init__(self, name=None):
        "Initialize an analytical YSO model"

        self.star = Star()
        self.disks = []
        self.envelopes = []

        self.accretion = False

        self.ambient = None

        self.minimum_temperature = None

        Model.__init__(self, name=name)

    def __setattr__(self, attribute, value):
        if attribute == 'accretion':
            if value is True:
                self.star.sources['uv'] = SphericalSource(name='uv', radius=self.star.radius)
                self.star.sources['xray'] = SphericalSource(name='xray', radius=self.star.radius)
            else:
                if 'uv' in self.star.sources:
                    self.star.sources.pop('uv')
                if 'xray' in self.star.sources:
                    self.star.sources.pop('xray')
        Model.__setattr__(self, attribute, value)

    # DENSITY COMPONENTS

    def add_ambient_medium(self):
        "Add an ambient medium in which the model is embedded"
        if self.ambient is not None:
            raise Exception("Ambient medium already present")
        ambient = AmbientMedium()
        self.ambient = ambient
        return ambient

    def add_flared_disk(self):
        "Add a flared disk to the geometry"
        disk = FlaredDisk()
        self.disks.append(disk)
        return disk

    def add_alpha_disk(self, definition='whitney'):
        "Add an alpha disk to the geometry"
        if definition == 'whitney':
            disk = AlphaDiskWhitney()
            disk.star = self.star
        else:
            raise Exception("Unknown alpha disk definition: %s" % definition)
        self.disks.append(disk)
        return disk

    def add_settled_disks(self, reference_disk, reference_size, eta=0.,
                          sizes=[], dust_files=[]):
        "Automatically create disks with varying degrees of settling"

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
            disk.h_0 *= (size / reference_size) ** (-eta)
            disk.dust = dust_files[i]
            self.disks.append(disk)

    def add_ulrich_envelope(self):
        "Add an infalling Ulrich envelope to the geometry"
        envelope = UlrichEnvelope()
        envelope.star = self.star
        self.envelopes.append(envelope)
        return envelope

    def add_power_law_envelope(self):
        "Add a power-law envelope to the geometry"
        envelope = PowerLawEnvelope()
        self.envelopes.append(envelope)
        return envelope

    def _check_all_set(self):
        for disk in self.disks:
            disk._check_all_set()
        for envelope in self.envelopes:
            envelope._check_all_set()
        if self.ambient is not None:
            self.ambient._check_all_set()

    # MIDPLANE OPTICAL DEPTH

    def print_midplane_tau(self, wavelength):
        for i, disk in enumerate(self.disks):
            if disk.mass > 0.:
                tau = disk.midplane_cumulative_density(np.array([disk.rmax])) \
                    * disk.dust.interp_chi_wav(wavelength)
                print "Disk %i: %.5e" % (i + 1, tau)

    def get_midplane_tau(self, r):

        self._check_all_set()

        # Find the total combined optical depth through the midplane
        tau_midplane = np.zeros(r.shape)

        # IMPLEMENT: PAHs

        for i, disk in enumerate(self.disks):
            if disk.mass > 0.:
                if disk.dust is None:
                    raise Exception("Disk %i dust not set" % i)
                nu, fnu = self.star.total_spectrum()
                tau_midplane += disk.midplane_cumulative_density(r) \
                              * disk.dust.optical_properties.chi_planck_spectrum(nu, fnu)

        for i, envelope in enumerate(self.envelopes):
            if envelope.exists():
                if envelope.dust is None:
                    raise Exception("envelope %i dust not set" % i)
                nu, fnu = self.star.total_spectrum()
                tau_midplane += envelope.midplane_cumulative_density(r) \
                              * envelope.dust.optical_properties.chi_planck_spectrum(nu, fnu)

        return tau_midplane

    def plot_midplane_tau(self, filename):

        tau_midplane = self.get_midplane_tau(self.grid.r_wall[1:])

        fig = mpl.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.loglog(self.grid.r[1:] / self.grid.r[1] - 1.,
                  tau_midplane[1:] - tau_midplane[:-1],
                  drawstyle='steps-mid')
        fig.savefig(filename)

    # COORDINATE GRID

    def radial_range(self):

        if len(self.disks) == 0 and len(self.envelopes) == 0:
            rmin = self.star.radius
        else:
            rmin_values = [disk.rmin for disk in self.disks] \
                        + [envelope.rmin for envelope in self.envelopes]
            if self.ambient is not None:
                rmin_values += [self.ambient.rmin]
            rmin = _min_none(*rmin_values)

        rmax_values = [self.star.radius]
        rmax_values += [disk.rmax for disk in self.disks] \
                     + [envelope.rmax for envelope in self.envelopes]
        if self.ambient is not None:
            rmax_values += [self.ambient.rmax]
        rmax = _max_none(*rmax_values)

        return rmin, rmax

    def set_cylindrical_polar_grid_auto(self, n_w, n_z, n_phi, zmax=None):
        self._set_polar_grid_auto(n_w, n_z, n_phi, 'cylindrical', zmax=zmax)

    def set_spherical_polar_grid_auto(self, n_r, n_theta, n_phi, rmax=None):
        self._set_polar_grid_auto(n_r, n_theta, n_phi, 'spherical', rmax=rmax)

    def _set_polar_grid_auto(self, n1, n2, n3, grid_type, zmax=None,
                             rmax=None):

        self.star._finalize()
        self._resolve_optically_thin_radii()

        if self.star.radius is None:
            raise Exception("The central source radius need to be defined before the grid can be set up")

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
            rmin_values = [disk.rmin for disk in self.disks] \
                        + [envelope.rmin for envelope in self.envelopes]
            if self.ambient is not None:
                rmin_values += [self.ambient.rmin]
            rmin = _min_none(*rmin_values)

        if not rmax:
            rmax_values = [2. * self.star.radius]
            rmax_values += [disk.rmax for disk in self.disks] \
                         + [envelope.rmax for envelope in self.envelopes]
            if self.ambient is not None:
                rmax_values += [self.ambient.rmax]
            rmax = _max_none(*rmax_values)

        if rmax < rmin:
            logger.warn("Grid rmax < rmin, model with consist only of central star")
            rmin = self.star.radius
            rmax = 2. * self.star.radius

        if np.isnan(rmin):
            raise Exception("R_min is NaN")

        if np.isnan(rmax):
            raise Exception("R_max is NaN")

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

        # Create interpolating function to find radius from optical depth
        # Can't do log interpolation because first element is zero
        r_interp = interp1d(tau_midplane, r)

        # Find where the second wall after rmin would be if grid cells were
        # spaced equally
        r_next_real = rmin * ((rmax / rmin) ** (1. / n_r) - 1.)

        # Find where the second wall would be if we put it at tau=0.1
        if tau_midplane[-1] <= 0.1:
            r_next_tau = rmax - rmin
        else:
            r_next_tau = r_interp(0.1) - rmin

        # Pick the smallest
        rnext = min(r_next_real, r_next_tau)

        # Make sure rnext isn't too small
        if rmin * (1. + 1.e-12) > rnext + rmin:
            rnext = rmin * 1.e-12

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
                zmin = min(zmin, disk.h_0 * (rmin / disk.r_0) ** disk.beta)

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
            self.set_spherical_polar_grid(r_wall, t_wall, p_wall)
        else:
            self.set_cylindrical_polar_grid(r_wall, z_wall, p_wall)

    # ACCRETION

    def setup_magnetospheric_accretion(self, mdot, rtrunc, fspot, xwav_min=0.001, xwav_max=0.01):
        '''
        Set up the model for magnetospheric accretion

        Parameters
        ----------
        mdot : float
            The accretion rate onto the star in cgs
        rtrunc:
            The magnetospheric truncation radius of the disk in cgs
        fspot : float
            The spot coverage fraction

        Notes
        -----
        This method currently assumes that the luminosity is split up into
        that which goes into the dust disk, gas disk, and the stellar surface.
        However, at this time the gas emission is not implemented.

        This method should be called once the stellar parameters have been
        otherwise initialized, and the disk parameters have to be set. This
        method cannot be called once the grid has been set, since this routine
        potentially changes the luminosity of the central source, potentially
        changing the dust sublimation radius.

        Calling this method causes the stellar parameters to be finalized,
        i.e. once this method has been called, none of the attributes of
        Model.star can be further modified.
        '''

        # For convenience
        lstar = self.star.sources['star'].luminosity

        # Tell the model that we are including accretion
        if not self.accretion:
            self.accretion = True

        if self.star.mass is None:
            raise Exception("Stellar mass is not set")

        # Find the luminosity dissipated in the shock
        lshock = G * self.star.mass * mdot * (1 / self.star.radius - 1 / rtrunc)

        # Hot spot parameters
        fluxratio = 0.5 * lshock / lstar / fspot
        teff = (lstar / (4. * pi * self.star.radius ** 2 * sigma)) ** 0.25  # Kelvin
        tshock = teff * (1 + fluxratio) ** 0.25  # Kelvin

        # Set the hot spot source
        self.star.sources['uv'].luminosity = lshock / 2. + lstar * fspot
        self.star.sources['uv'].temperature = tshock

        # X-rays from 0.1 to 10nm
        wav = np.logspace(np.log10(xwav_min), np.log10(xwav_max), 100)[::-1]
        nu = c * 1.e4 / wav
        fnu = np.repeat(1., nu.shape)

        # Set the X-ray source
        self.star.sources['xray'].luminosity = lshock / 2.
        self.star.sources['xray'].spectrum = (nu, fnu)

        # Reduce the total luminosity from the original source
        self.star.sources['star'].luminosity *= 1 - fspot

    # RESOLVERS

    def _resolve_optically_thin_radii(self):
        if not self.star.isfinal():
            raise Exception("Stellar parameters need to be finalized before resolving radiation-dependent radii")
        for i, disk in enumerate(self.disks):
            if isinstance(disk.rmin, OptThinRadius):
                if disk.dust is None:
                    raise Exception("Disk %i dust not set" % (i + 1))
                disk.rmin = disk.rmin.evaluate(self.star, disk.dust)
            if isinstance(disk.rmax, OptThinRadius):
                if disk.dust is None:
                    raise Exception("Disk %i dust not set" % (i + 1))
                disk.rmax = disk.rmax.evaluate(self.star, disk.dust)
        for i, envelope in enumerate(self.envelopes):
            if isinstance(envelope.rmin, OptThinRadius):
                if envelope.dust is None:
                    raise Exception("Envelope %i dust not set" % (i + 1))
                envelope.rmin = envelope.rmin.evaluate(self.star, envelope.dust)
            if isinstance(envelope.rmax, OptThinRadius):
                if envelope.dust is None:
                    raise Exception("Envelope %i dust not set" % (i + 1))
                envelope.rmax = envelope.rmax.evaluate(self.star, envelope.dust)
        if self.ambient is not None:
            if isinstance(self.ambient.rmin, OptThinRadius):
                if self.ambient.dust is None:
                    raise Exception("Ambient medium dust not set")
                self.ambient.rmin = self.ambient.rmin.evaluate(self.star, self.ambient.dust)
            if isinstance(self.ambient.rmax, OptThinRadius):
                if self.ambient.dust is None:
                    raise Exception("Ambient medium dust not set")
                self.ambient.rmax = self.ambient.rmax.evaluate(self.star, self.ambient.dust)

    # OUTPUT

    def set_minimum_temperature(self, temperature):
        self.minimum_temperature = temperature

    def write(self, filename=None, compression=True, copy=True,
              absolute_paths=False, wall_dtype=float,
              physics_dtype=float, overwrite=True, merge_if_possible=True):
        '''
        Write the model input parameters to an HDF5 file

        Parameters
        ----------
        filename: str
            The name of the input file to write. If no name is specified, the
            filename is constructed from the model name.
        compression: bool
            Whether to compress the datasets inside the HDF5 file.
        copy: bool
            Whether to copy all external content into the input file, or
            whether to just link to external content.
        absolute_paths: bool
            If copy=False, then if absolute_paths is True, absolute filenames
            are used in the link, otherwise the path relative to the input
            file is used.
        wall_dtype: type
            Numerical type to use for wall positions.
        physics_dtype: type
            Numerical type to use for physical grids.
        overwrite: bool
            Whether to overwrite any pre-existing file
        merge_if_possible: bool
            Whether to merge densities for the same dust populations
        '''

        if self.grid is None:
            raise Exception("The coordinate grid needs to be defined before calling AnalyticalModelYSO.write(...)")

        self.grid.reset_quantities()

        for i, disk in enumerate(self.disks):

            if disk.rmin >= disk.rmax:
                logger.warn("Disk rmin >= rmax, ignoring density contribution")
            elif disk.mass == 0.:
                logger.warn("Disk mass is zero, ignoring density contribution")
            else:

                if not disk.dust:
                    raise Exception("Disk %i dust not set" % (i + 1))
                self.add_density_grid(disk.density(self.grid), disk.dust,
                                      minimum_temperature=self.minimum_temperature,
                                      merge_if_possible=merge_if_possible)

        for i, envelope in enumerate(self.envelopes):

            if envelope.rmin >= envelope.rmax:
                logger.warn("Envelope rmin >= rmax, ignoring density contribution")
            elif isinstance(envelope, UlrichEnvelope) and envelope.rho_0 == 0.:
                logger.warn("Ulrich envelope has zero density everywhere, ignoring density contribution")
            elif isinstance(envelope, PowerLawEnvelope) and envelope.mass == 0.:
                logger.warn("Power-law envelope has zero density everywhere, ignoring density contribution")
            else:

                if not envelope.dust:
                    raise Exception("Envelope dust not set")
                self.add_density_grid(envelope.density(self.grid), envelope.dust,
                                      minimum_temperature=self.minimum_temperature,
                                      merge_if_possible=merge_if_possible)

                if envelope.cavity is not None:
                    if envelope.cavity.theta_0 == 0.:
                        logger.warn("Cavity opening angle is zero, ignoring density contribution")
                    elif envelope.cavity.rho_0 == 0.:
                        logger.warn("Cavity density is zero, ignoring density contribution")
                    else:
                        if not envelope.cavity.dust:
                            raise Exception("Cavity dust not set")
                        self.add_density_grid(envelope.cavity.density(self.grid), envelope.cavity.dust,
                                              minimum_temperature=self.minimum_temperature,
                                              merge_if_possible=merge_if_possible)

        # AMBIENT MEDIUM

        if self.ambient is not None:

            if self.ambient.density == 0.:
                logger.warn("Ambient medium has zero density, ignoring density contribution")
            else:

                ambient = self.ambient

                if not ambient.dust:
                    raise Exception("Ambient medium dust not set")

                # Find the density of the ambient medium
                density_amb = ambient.density(self.grid)

                if len(self.density) > 0:

                    # Find total density in other components
                    shape = list(self.grid.shape)
                    shape.insert(0, len(self.density))
                    density_sum = np.sum(np.vstack(self.density).reshape(*shape), axis=0)

                    density_amb -= density_sum
                    density_amb[density_amb < 0.] = 0.

                self.add_density_grid(density_amb, ambient.dust,
                                      minimum_temperature=self.minimum_temperature,
                                      merge_if_possible=merge_if_possible)

        # SOURCES

        # Star

        if self.star.sources['star'].luminosity > 0:
            if self.star.sources['star'] not in self.sources:
                self.add_source(self.star.sources['star'])

        # Accretion

        if self.accretion:

            if self.star.sources['uv'].luminosity > 0.:
                if self.star.sources['uv'] not in self.sources:
                    self.add_source(self.star.sources['uv'])

            if self.star.sources['xray'].luminosity > 0.:
                if self.star.sources['xray'] not in self.sources:
                    self.add_source(self.star.sources['xray'])

            for i, disk in enumerate(self.disks):

                if disk.rmin >= disk.rmax:
                    logger.warn("Disk rmin >= rmax, ignoring accretion luminosity")
                elif disk.mass == 0.:
                    logger.warn("Disk mass is zero, ignoring accretion luminosity")
                elif disk.lvisc == 0.:
                    logger.warn("Disk viscous luminosity is zero, ignoring accretion luminosity")
                else:
                    self.add_map_source(luminosity=disk.lvisc, map=disk.accretion_luminosity(self.grid), name='accdisk%i' % i)

        Model.write(self, filename=filename, compression=compression,
                    copy=copy, absolute_paths=absolute_paths,
                    wall_dtype=wall_dtype, physics_dtype=physics_dtype,
                    overwrite=overwrite)
