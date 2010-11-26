import warnings
from copy import deepcopy

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as mpl

from hyperion.model import Model
from hyperion.densities import FlaredDisk, PowerLawEnvelope, UlrichEnvelope, AmbientMedium
from hyperion.util.interpolate import interp1d_fast_loglog
from hyperion.util.constants import pi, sigma, c
from hyperion.sources import SphericalSource, SpotSource
from hyperion.util.functions import FreezableClass
from hyperion.dust import SphericalDust
from hyperion.util.convenience import OptThinRadius


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
        self._freeze()

    def add_spot(self, *args, **kwargs):
        self.sources['star'].add_spot(SpotSource(*args, **kwargs))

    def __setattr__(self, attribute, value):
        if self.isfinal():
            raise Exception("Attribute %s can no longer be changed" % attribute)
        if attribute in ['luminosity', 'temperature', 'spectrum']:
            self.sources['star'].__setattr__(attribute, value)
            return
        elif attribute == 'radius':
            for source in self.sources:
                self.sources[source].radius = value
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

    def __init__(self, name):
        "Initialize an analytical YSO model"

        self.star = Star()
        self.disks = []
        self.envelopes = []

        self.accretion = False

        self.ambient = None

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
        disk.star = self.star
        self.disks.append(disk)
        return disk

    def add_settled_disks(self, reference_disk, reference_size, eta=0.,
                          sizes=[], dust_files=[]):
        "Automatically create disks with varying degrees of settling"

        exists = False

        for disk in self.disks:
            if disk is reference_disk:
                warnings.warn("Reference disk already exists, not re-adding")
                exists = True

        if not exists:
            warnings.warn("Reference disk does not exist, adding")
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
        envelope.star = self.star
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
                dust = SphericalDust(disk.dust)
                tau = disk.midplane_cumulative_density(np.array([disk.rmax])) \
                    * dust.interp_chi_wav(wavelength)
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
                dust = SphericalDust(disk.dust)
                nu, fnu = self.star.total_spectrum()
                tau_midplane += disk.midplane_cumulative_density(r) \
                              * dust.chi_planck_spectrum(nu, fnu)

        for i, envelope in enumerate(self.envelopes):
            if envelope.exists():
                if envelope.dust is None:
                    raise Exception("envelope %i dust not set" % i)
                dust = SphericalDust(envelope.dust)
                nu, fnu = self.star.total_spectrum()
                tau_midplane += envelope.midplane_cumulative_density(r) \
                              * dust.chi_planck_spectrum(nu, fnu)

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

    def set_cylindrical_polar_grid_auto(self, n_w, n_z, n_phi, zmax=None):
        self._set_polar_grid_auto(n_w, n_z, n_phi, 'cylindrical', zmax=zmax)

    def set_spherical_polar_grid_auto(self, n_r, n_theta, n_phi, rmax=None):
        self._set_polar_grid_auto(n_r, n_theta, n_phi, 'spherical', rmax=rmax)

    def _set_polar_grid_auto(self, n1, n2, n3, grid_type, zmax=None,
                             rmax=None):

        self.star._finalize()
        self._resolve_optically_thin_radii()

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
            warnings.warn("Grid rmax < rmin, model with consist only of central star")
            rmin = self.star.radius
            rmax = 2. * self.star.radius

        # RADIAL WALLS

        # Set first wall to be at half the stellar radius to avoid any
        # numerical problems. Set second wall to be the disk inner radius.

        r_wall = [self.star.radius / 2., rmin]

        # Define a radial grid to compute the midplane column density on
        r = np.logspace(-20., np.log10(rmax / rmin), 100000) * rmin + rmin

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
            r_next_tau = rmax
        else:
            r_next_tau = r_interp(0.1) - rmin

        # Pick the smallest
        rnext = min(r_next_real, r_next_tau)

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

    def setup_magnetospheric_accretion(self, lacc, rtrunc, fspot, disk):
        '''
        Set up the model for magnetospheric accretion

        Parameters
        ----------
        lacc : float
            The total luminosity to be released via accretion
        rtrunc : float
            The magnetospheric truncation radius
        fspot : float
            The spot coverage fraction
        disk : FlaredDisk
            The disk that will emit the viscously dissipated accretion energy

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
        self.accretion = True

        # Find the total luminosity of accretion shock on star
        frac_shock = 1 / self.star.radius - 1 / rtrunc
        frac_disk = 0.5 / rtrunc - 0.5 / disk.rmax

        # Find the luminosity dissipated in the shock
        lshock = lacc * frac_shock / (frac_shock + frac_disk)

        # Hot spot parameters
        fluxratio = 0.5 * lshock / lstar / fspot
        teff = (lstar / (4. * pi * self.star.radius ** 2 * sigma)) ** 0.25  # Kelvin
        tshock = teff * (1 + fluxratio) ** 0.25  # Kelvin

        # Set the hot spot source
        self.star.sources['uv'].luminosity = lshock / 2. + lstar * fspot
        self.star.sources['uv'].temperature = tshock

        # X-rays from 0.1 to 10nm
        wav = np.logspace(-3., -2., 100)[::-1]
        nu = c * 1.e4 / wav
        fnu = np.repeat(1., nu.shape)

        # Set the X-ray source
        self.star.sources['xray'].luminosity = lshock / 2.
        self.star.sources['xray'].spectrum = (nu, fnu)

        # Reduce the total luminosity from the original source
        self.star.sources['star'].luminosity *= 1 - fspot

        # Prevent any further modifiations to the star
        self.star._finalize()

        # Resolve any sublimation radii in envelopes/disks

        # Set luminosity from viscous dissipation in disk
        # For this, only find the part that is inside the dust disk
        disk.lvisc = lacc  # WRONG!!

    # RESOLVERS

    def _resolve_optically_thin_radii(self):
        if not self.star.isfinal():
            raise Exception("Stellar parameters need to be finalized before resolving radiation-dependent radii")
        for disk in self.disks:
            if isinstance(disk.rmin, OptThinRadius):
                disk.rmin = disk.rmin.evaluate(self.star, SphericalDust(disk.dust))
            if isinstance(disk.rmax, OptThinRadius):
                disk.rmax = disk.rmax.evaluate(self.star, SphericalDust(disk.dust))
        for envelope in self.envelopes:
            if isinstance(envelope.rmin, OptThinRadius):
                envelope.rmin = envelope.rmin.evaluate(self.star, SphericalDust(envelope.dust))
            if isinstance(envelope.rmax, OptThinRadius):
                envelope.rmax = envelope.rmax.evaluate(self.star, SphericalDust(envelope.dust))
        if self.ambient is not None:
            if isinstance(self.ambient.rmin, OptThinRadius):
                self.ambient.rmin = self.ambient.rmin.evaluate(self.star, SphericalDust(self.ambient.dust))
            if isinstance(self.ambient.rmax, OptThinRadius):
                self.ambient.rmax = self.ambient.rmax.evaluate(self.star, SphericalDust(self.ambient.dust))

    # OUTPUT

    def write(self, **kwargs):

        self.reset_density()
        self.reset_sources()

        for i, disk in enumerate(self.disks):

            if disk.mass > 0.:
                if not disk.dust:
                    raise Exception("Disk %i dust not set" % (i + 1))
                self.add_density_grid(disk.density(self.grid), disk.dust)

        for i, envelope in enumerate(self.envelopes):

            if envelope.exists():
                if not envelope.dust:
                    raise Exception("Envelope dust not set")
                self.add_density_grid(envelope.density(self.grid),
                                      envelope.dust)

            if envelope.cavity is not None:
                if envelope.cavity.exists():
                    if not envelope.cavity.dust:
                        raise Exception("Cavity dust not set")
                    self.add_density_grid(envelope.cavity.density(self.grid),
                                          envelope.cavity.dust)

        # AMBIENT MEDIUM

        if self.ambient is not None:

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

            self.add_density_grid(density_amb, ambient.dust)

        # SOURCES

        # Star

        if self.star.sources['star'].luminosity > 0:
            self.add_source(self.star.sources['star'])

        # Accretion

        if self.accretion:

            if self.star.sources['uv'].luminosity > 0.:
                self.add_source(self.star.sources['uv'])

            if self.star.sources['xray'].luminosity > 0.:
                self.add_source(self.star.sources['xray'])

            for i, disk in enumerate(self.disks):
                self.add_map_source(luminosity=disk.lvisc, map=disk.accretion_luminosity(self.grid), name='accdisk%i' % i)

        Model.write(self, **kwargs)
