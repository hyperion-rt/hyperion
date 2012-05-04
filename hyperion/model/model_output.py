from __future__ import print_function, division

import os

import h5py
import numpy as np

from ..util.constants import c, pi
from ..util.functions import FreezableClass
from ..dust import SphericalDust
from ..util.logger import logger
from ..util.decorator import decorator

from .helpers import find_last_iteration

STOKESD = {}
STOKESD['I'] = 0
STOKESD['Q'] = 1
STOKESD['U'] = 2
STOKESD['V'] = 3

LABEL = {}
LABEL['I'] = '$\lambda\, F_\lambda$'
LABEL['Q'] = '$\lambda\, F_\lambda$ [stokes=Q]'
LABEL['U'] = '$\lambda\, F_\lambda$ [stokes=U]'
LABEL['V'] = '$\lambda\, F_\lambda$ [stokes=V]'
LABEL['linpol'] = "Total linear polarization fraction"
LABEL['circpol'] = "Total circular polarization fraction"

UNITS_LABEL = {}
UNITS_LABEL['ergs/s'] = '(ergs/s)'
UNITS_LABEL['ergs/cm^2/s'] = '(ergs/cm$^2$/s)'
UNITS_LABEL['ergs/cm^2/s/Hz'] = '(ergs/cm$^2$/s/Hz)'
UNITS_LABEL['Jy'] = 'Jy'
UNITS_LABEL['mJy'] = 'mJy'
UNITS_LABEL['MJy/sr'] = 'MJy/sr'


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

# We now define a decorator for methods that needs access to the output HDF5
# file. This is necessary because h5py has issues with links pointing to
# groups that are in open files.


def on_the_fly_hdf5(f):
    return decorator(_on_the_fly_hdf5, f)


def _on_the_fly_hdf5(f, *args, **kwargs):
    preset = args[0].file is not None
    if not preset:
        args[0].file = h5py.File(args[0].filename, 'r')
    results = f(*args, **kwargs)
    if not preset:
        args[0].file.close()
        args[0].file = None
    return results


class ModelOutput(FreezableClass):

    def __init__(self, filename):
        '''
        Create a ModelOutput instance that can be used to access the
        data in the output file.

        Parameters
        ----------
        name : str
            The name of the model output file.
        '''

        # Check that file exists
        if not os.path.exists(filename):
            raise IOError("File not found: %s" % filename)

        # Open file and store handle to object
        # (but don't read in the contents yet)
        self.filename = filename
        self.file = None

    @on_the_fly_hdf5
    def get_sed(self, stokes='I', group=0, technique='peeled',
                distance=None, component='total', inclination='all',
                aperture='all', uncertainties=False, units=None,
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
            The peeloff group (zero-based). If multiple peeloff image groups
            were requested, this can be used to select between them. The
            default is to return the first group. This option is only used if
            technique='peeled'.

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
            The default is ``'ergs/cm^2/s'``. If a distance is not specified,
            then this option is ignored, and the output units are ergs/s.

        source_id, dust_id : int or str, optional
            If the output file was made with track_origin='detailed', a
            specific source and dust component can be specified. If 'all' is
            specified, then all components are returned individually. If
            neither of these are not specified, then the total component
            requested for all sources or dust types is returned.

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

        # Set default units
        if units is None:
            if distance is None:
                units = 'ergs/s'
            else:
                units = 'ergs/cm^2/s'

        # Check argument types
        if type(stokes) is not str:
            raise ValueError("stokes argument should be a string")

        # Check for inconsistent parameters
        if distance is not None and stokes in ['linpol', 'circpol']:
            raise Exception("Cannot scale linear or circular polarization degree by distance")

        if technique == 'peeled':
            n_groups = len(self.file['Peeled'])
            if (group < 0 and -group <= n_groups) or (group >= 0 and group < n_groups):
                g = self.file['Peeled/group_%05i' % (group + 1)]
            else:
                raise ValueError('File only contains %i image/SED group(s)' % n_groups)
        else:
            g = self.file['Binned']

        if not 'seds' in g:
            raise Exception("Group %i does not contain any SEDs" % group)

        # Check that uncertainties are present if requested
        if uncertainties and not 'seds_unc' in g:
            raise Exception("Uncertainties requested but not present in file")

        if 'track_origin' in g['seds'].attrs:

            track_origin = g['seds'].attrs['track_origin'].decode('utf-8')

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

                    io = ((io - (io + 1) % 2 + 1) * ns + (io - io % 2) * nd) / 2

                    if component.startswith('source'):
                        if source_id is None or source_id == 'all':
                            io = (io, io + ns)
                        else:
                            if source_id < 1 or source_id > ns:
                                raise Exception("source_id should be between 1 and %i" % ns)
                            io = io + (source_id - 1)
                    else:
                        if dust_id is None or dust_id == 'all':
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
            wav = np.logspace(np.log10(wavmax), np.log10(wavmin), g['seds'].shape[-1] * 2 + 1)[1::2]
            nu = c / wav * 1.e4
        else:
            nu = g['frequencies']['nu']
            wav = c / nu * 1.e4

        flux = g['seds'].value
        if uncertainties:
            unc = g['seds_unc'].value

        try:
            inside_observer = g.attrs['inside_observer'].decode('utf-8').lower() == 'yes'
        except:
            inside_observer = False

        if inside_observer and distance is not None:
            raise ValueError("Cannot specify distance for inside observers")

        # Optionally scale by distance
        if distance is not None or inside_observer:

            # Convert to the correct units
            if units == 'ergs/cm^2/s':
                scale = np.repeat(1., len(nu))
            elif units == 'ergs/cm^2/s/Hz':
                scale = 1. / nu
            elif units == 'Jy':
                scale = 1.e23 / nu
            elif units == 'mJy':
                scale = 1.e26 / nu
            else:
                raise Exception("Unknown units: %s" % units)

            # Scale by distance
            if distance:
                scale *= 1. / (4. * pi * distance ** 2)

        else:

            if units != 'ergs/s':
                raise ValueError("Since distance= is not specified, units should be set to ergs/s")

            # Units here are not technically ergs/cm^2/s but ergs/s
            scale = np.repeat(1., len(nu))

        # If in 32-bit mode, need to convert to 64-bit because of scaling/polarization to be safe
        if flux.dtype == np.float32:
            flux = flux.astype(np.float64)
        if uncertainties and unc.dtype == np.float32:
            unc = unc.astype(np.float64)

        # If a stokes component is requested, scale the images. Frequency is
        # the last dimension, so this compact notation can be used.

        if stokes in STOKESD:
            flux *= scale
            if uncertainties:
                unc *= scale

        # We now slice the SED array to end up with what the user requested.
        # Note that we slice from the last to the first dimension to ensure that
        # we always know what the slices are. In this section, we make use of
        # the fact that when writing array[i] with a multi-dimensional array,
        # the index i applies only to the first dimension. So flux[1] really
        # means flux[1, :, :, :, :].

        if aperture == 'all':
            pass
        else:
            flux = flux[:, :, :, aperture]
            if uncertainties:
                unc = unc[:, :, :, aperture]

        if inclination == 'all':
            pass
        else:
            flux = flux[:, :, inclination]
            if uncertainties:
                unc = unc[:, :, inclination]

        # Select correct origin component

        if component == 'total':
            flux = np.sum(flux, axis=1)
            if uncertainties:
                unc = np.sqrt(np.sum(unc ** 2, axis=1))
        elif component in ['source_emit', 'dust_emit', 'source_scat', 'dust_scat']:
            if type(io) is tuple:
                start, end = io
                flux = flux[:, start:end]
                if uncertainties:
                    unc = unc[:, start:end]
                if (component.startswith('source') and source_id is None) or \
                   (component.startswith('dust') and dust_id is None):
                    flux = np.sum(flux, axis=1)
                    if uncertainties:
                        unc = np.sqrt(np.sum(unc ** 2, axis=1))
            else:
                flux = flux[:, io]
                if uncertainties:
                    unc = unc[:, io]
        else:
            raise Exception("Unknown component: %s" % component)

        # Select correct Stokes component

        if stokes in STOKESD:
            flux = flux[STOKESD[stokes]]
            if uncertainties:
                unc = unc[STOKESD[stokes]]
        elif stokes == 'linpol':
            if uncertainties:
                f = np.vectorize(mc_linear_polarization)
                flux, unc = f(flux[0], unc[0], flux[1], unc[1], flux[2], unc[2])
            else:
                flux = np.sqrt((flux[1] ** 2 + flux[2] ** 2) / flux[0] ** 2)
                flux[np.isnan(flux)] = 0.
        elif stokes == 'circpol':
            if uncertainties:
                f = np.vectorize(mc_circular_polarization)
                flux, unc = f(flux[0], unc[0], flux[3], unc[3])
            else:
                flux = np.abs(flux[3] / flux[0])
                flux[np.isnan(flux)] = 0.
        else:
            raise ValueError("Unknown Stokes parameter: %s" % stokes)

        if uncertainties:
            return wav, flux, unc
        else:
            return wav, flux


    def plot_sed(self, axes=None, filename=None,
                 wmin=0.01, wmax=5000., fmin=None, fmax=None,
                 color='black', labels=True, **kwargs):
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

        wmin, wmax : float, optional
            The range in wavelengths to show (in microns).

        fmin, fmax : float, optional
            The range in fluxes to show (in ergs/s if distance is not
            specified, or in ergs/cm^2/s if distance is specified).

        color : str, optional
            The color to plot the SED(s) in.

        labels : bool, optional
            Whether or not to show the axis labels

        All additional parameters are passed to get_sed, and any remaining
        parameters are passed to Axes.plot, Axes.loglog, or Axes.errorbar
        (depending on which one is used)

        Returns
        -------

        axes : matplotlib.pyplot.Axes instance (if `axes` was set)
            The updated matplotlib Axes

        figure : matplotlib.pyplot.Figure instance (if `axes` and `filename` were not set)
            The matplotlib Figure created

        '''

        import matplotlib.pyplot as plt

        # TODO: make this routine work when origin='detailed' and source_id or
        # dust_id are 'all'

        # Check for inconsistent parameters
        if axes is not None and filename is not None:
            raise Exception("Cannot specify both an axes instance and a filename")

        if axes is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = axes

        # Store a few of the kwargs into local variables
        if 'inclination' in kwargs:
            inclination = kwargs['inclination']
        else:
            inclination = 'all'

        if 'aperture' in kwargs:
            aperture = kwargs['aperture']
        else:
            aperture = 'all'

        if 'uncertainties' in kwargs:
            uncertainties = kwargs['uncertainties']
        else:
            uncertainties = False

        if 'stokes' not in kwargs:
            kwargs['stokes'] = 'I'

        # Make sure that SED cube always has three dimensions, so need to request inclinations and apertures in a list
        if type(inclination) is not str and np.isscalar(inclination):
            inclination = [inclination]
        if type(aperture) is not str and np.isscalar(aperture):
            aperture = [aperture]

        seds = self.get_sed(**kwargs)

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
                    logger.warn("All fluxes are zero - cannot plot in log-log plot")
                else:

                    # Plot main SED
                    if kwargs['stokes'] in ['Q', 'U', 'V']:
                        ax.plot(wav, nufnu[ii, ia, :], color=color)
                        ax.set_xscale('log')
                    else:
                        ax.loglog(wav, nufnu[ii, ia, :], color=color)

                    # Check whether the maximum flux needs to be changed
                    maxflux = max(maxflux, np.nanmax(nufnu))

                    if uncertainties:
                        nufnu_lower = np.maximum(nufnu[ii, ia, :] * 1.e-10, nufnu[ii, ia, :] - nufnu_unc[ii, ia, :])
                        nufnu_unc_lower = nufnu[ii, ia, :] - nufnu_lower
                        ax.errorbar(wav, nufnu[ii, ia, :], yerr=[nufnu_unc_lower, nufnu_unc[ii, ia, :]], fmt=None, ecolor=color, **kwargs)

        if not fmax:
            if kwargs['stokes'] in ['Q', 'U', 'V']:
                fmax = maxflux * 1.5
            elif kwargs['stokes'] in ['linpol', 'circpol']:
                fmax = 1.
            else:
                fmax = maxflux * 10.

        if not fmin:
            if kwargs['stokes'] in ['Q', 'U', 'V']:
                fmin = -fmax
            else:
                fmin = fmax * 1.e-6

        ax.set_xlim(wmin, wmax)
        ax.set_ylim(fmin, fmax)

        if labels:
            ax.set_xlabel('$\lambda$ ($\mu$m)')
            if kwargs['stokes'] in ['I', 'Q', 'U', 'V']:
                if kwargs['units'] is None:
                    if kwargs['distance'] is None:
                        ax.set_ylabel(LABEL[kwargs['stokes']] + " " + UNITS_LABEL['ergs/s'])
                    else:
                        ax.set_ylabel(LABEL[kwargs['stokes']] + " " + UNITS_LABEL['ergs/cm^2/s'])
                else:
                    ax.set_ylabel(LABEL[kwargs['stokes']] + " " + UNITS_LABEL[kwargs['units']])
            else:
                ax.set_ylabel(LABEL[kwargs['stokes']])

        if filename:
            fig.savefig(filename)
            plt.close(fig)
            return
        elif axes is None:
            return fig
        else:
            return ax


    @on_the_fly_hdf5
    def get_image(self, stokes='I', group=0, technique='peeled',
                  distance=None, component='total', inclination='all',
                  uncertainties=False, units=None,
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
            The peeloff group (zero-based). If multiple peeloff image groups
            were requested, this can be used to select between them. The
            default is to return the first group. This option is only used if
            technique='peeled'.

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

        source_id, dust_id : int or str, optional
            If the output file was made with track_origin='detailed', a
            specific source and dust component can be specified. If 'all' is
            specified, then all components are returned individually. If
            neither of these are not specified, then the total component
            requested for all sources or dust types is returned.

        All additional parameters are passed to get_image.

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

        # Set default units
        if units is None:
            if distance is None:
                units = 'ergs/s'
            else:
                units = 'ergs/cm^2/s'

        # Check argument types
        if type(stokes) is not str:
            raise ValueError("stokes argument should be a string")

        # Check for inconsistent parameters
        if distance is not None and stokes in ['linpol', 'circpol']:
            raise Exception("Cannot scale linear or circular polarization degree by distance")

        if technique == 'peeled':
            n_groups = len(self.file['Peeled'])
            if (group < 0 and -group <= n_groups) or (group >= 0 and group < n_groups):
                g = self.file['Peeled/group_%05i' % (group + 1)]
            else:
                raise ValueError('File only contains %i image/SED group(s)' % n_groups)
        else:
            g = self.file['Binned']

        if not 'images' in g:
            raise Exception("Group %i does not contain any images" % group)

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
                        if source_id is None or source_id == 'all':
                            io = (io, io + ns)
                        else:
                            if source_id < 1 or source_id > ns:
                                raise Exception("source_id should be between 1 and %i" % ns)
                            io = io + (source_id - 1)
                    else:
                        if dust_id is None or dust_id == 'all':
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
            wav = np.logspace(np.log10(wavmax), np.log10(wavmin), g['images'].shape[-1] * 2 + 1)[1::2]
            nu = c / wav * 1.e4
        else:
            nu = g['frequencies']['nu']
            wav = c / nu * 1.e4

        flux = g['images'].value
        if uncertainties:
            unc = g['images_unc'].value

        try:
            inside_observer = g.attrs['inside_observer'].decode('utf-8').lower() == 'yes'
        except:
            inside_observer = False

        if inside_observer and distance is not None:
            raise ValueError("Cannot specify distance for inside observers")

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
                ny, nx = flux.shape[-2:]

                # Find pixel resolution in radians/pixel
                if inside_observer:
                    pix_dx = abs(np.radians(xmax - xmin) / float(nx))
                    pix_dy = abs(np.radians(ymax - ymin) / float(ny))
                else:
                    pix_dx = abs(np.arctan((xmax - xmin) / float(nx) / distance))
                    pix_dy = abs(np.arctan((ymax - ymin) / float(ny) / distance))

                # Find pixel area in steradians
                pix_area = pix_dx * pix_dy

                scale = 1.e17 / nu / pix_area

            else:
                raise Exception("Unknown units: %s" % units)

            # Scale by distance
            if distance:
                scale *= 1. / (4. * pi * distance ** 2)

        else:

            if units != 'ergs/s':
                raise ValueError("Since distance= is not specified, units should be set to ergs/s")

            scale = np.repeat(1., len(nu))

        # If in 32-bit mode, need to convert to 64-bit because of scaling/polarization to be safe
        if flux.dtype == np.float32:
            flux = flux.astype(np.float64)
        if uncertainties and unc.dtype == np.float32:
            unc = unc.astype(np.float64)

        # If a stokes component is requested, scale the images. Frequency is
        # the last dimension, so this compact notation can be used.

        if stokes in STOKESD:
            flux *= scale
            if uncertainties:
                unc *= scale

        # We now slice the image array to end up with what the user requested.
        # Note that we slice from the last to the first dimension to ensure that
        # we always know what the slices are. In this section, we make use of
        # the fact that when writing array[i] with a multi-dimensional array,
        # the index i applies only to the first dimension. So flux[1] really
        # means flux[1, :, :, :, :].

        if inclination == 'all':
            pass
        else:
            flux = flux[:, :, inclination]
            if uncertainties:
                unc = unc[:, :, inclination]

        # Select correct origin component

        if component == 'total':
            flux = np.sum(flux, axis=1)
            if uncertainties:
                unc = np.sqrt(np.sum(unc ** 2, axis=1))
        elif component in ['source_emit', 'dust_emit', 'source_scat', 'dust_scat']:
            if type(io) is tuple:
                start, end = io
                flux = flux[:, start:end]
                if uncertainties:
                    unc = unc[:, start:end]
                if (component.startswith('source') and source_id is None) or \
                   (component.startswith('dust') and dust_id is None):
                    flux = np.sum(flux, axis=1)
                    if uncertainties:
                        unc = np.sqrt(np.sum(unc ** 2, axis=1))
            else:
                flux = flux[:, io]
                if uncertainties:
                    unc = unc[:, io]
        else:
            raise Exception("Unknown component: %s" % component)

        # Select correct Stokes component

        if stokes in STOKESD:
            flux = flux[STOKESD[stokes]]
            if uncertainties:
                unc = unc[STOKESD[stokes]]
        elif stokes == 'linpol':
            if uncertainties:
                f = np.vectorize(mc_linear_polarization)
                flux, unc = f(flux[0], unc[0], flux[1], unc[1], flux[2], unc[2])
            else:
                flux = np.sqrt((flux[1] ** 2 + flux[2] ** 2) / flux[0] ** 2)
                flux[np.isnan(flux)] = 0.
        elif stokes == 'circpol':
            if uncertainties:
                f = np.vectorize(mc_circular_polarization)
                flux, unc = f(flux[0], unc[0], flux[3], unc[3])
            else:
                flux = np.abs(flux[3] / flux[0])
                flux[np.isnan(flux)] = 0.
        else:
            raise ValueError("Unknown Stokes parameter: %s" % stokes)

        if uncertainties:
            return wav, flux, unc
        else:
            return wav, flux


    def plot_image(self, wavelength, axes=None, filename=None, vmin=None,
                   vmax=None, cmap=None, labels=True, **kwargs):
        '''
        Plot an image

        Parameters
        ----------

        wavelength : float
            The wavelength to plot the image for, in microns. The
            wavelength closest to this will be chose.

        axes : matplotlib.pyplot.Axes instance, optional
            The matplotlib Axes to plot the SED(s) into (incompatible with
            the `filename` option).

        filename : str, optional
            The file to write the plot to (incompatible with the `ax`
            option).

        vmin, vmax : float, optional
            The range in fluxes to show (in ergs/s if distance is not
            specified, or in ergs/cm^2/s if distance is specified).

        cmap : matplotlib colormap
            The color to plot the SED(s) in.

        labels : bool, optional
            Whether or not to show the axis labels

        All additional parameters are passed to get_image.

        Returns
        -------

        axes : matplotlib.pyplot.Axes instance (if `axes` was set)
            The updated matplotlib Axes are returned if they were passed in
            as input.

        figure : matplotlib.pyplot.Figure instance (if `axes` and `filename` were not set)
            The matplotlib Figure created
        '''

        import matplotlib.pyplot as plt

        if axes is None:
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        else:
            ax = axes

        if 'uncertainties' in kwargs and kwargs['uncertainties']:
            kwargs['uncertainties'] = False
            logger.warn("plot_sed does not support uncertainties at this time, setting uncertainties=False")

        # Retrieved the necessary image
        wav, nufnu = self.get_image(**kwargs)

        # Find index closest to wavelength requested
        iw = np.argmin(np.abs(wav - wavelength))

        # Plot the image
        ax.imshow(nufnu[iw, :, :], vmin=vmin, vmax=vmax, cmap=cmap,
                  aspect='auto')

        # Save figure if necessary
        if filename:
            fig.savefig(filename)
            plt.close(fig)
            return
        elif axes is None:
            return fig
        else:
            return ax


    @on_the_fly_hdf5
    def get_available_components(self, iteration=-1):
        '''
        Find out what physical components are available in the output file

        Parameters
        ----------

        iteration : integer, optional
            The iteration to retrieve the grid for. The default is to return the components for the last iteration
        '''

        # If iteration is last one, find iteration number
        if iteration == -1:
            iteration = find_last_iteration(self.file)

        # Return components
        components = list(self.file['iteration_%05i' % iteration].keys())
        if 'specific_energy' in components:
            components.append('temperature')
        return components


    @on_the_fly_hdf5
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

        dust_id : 'all' or int
            If applicable, the ID of the dust type to extract the grid for (does not apply to n_photons)


        Returns
        -------

        array : numpy.array instance for regular grids
            The physical grid in cgs.

        Notes
        -----

        At the moment, this method only works on regular grids, not AMR or Oct-tree grids
        '''

        # Check that dust_id was not specified if grid is n_photons
        if name == 'n_photons':
            if dust_id != 'all':
                raise ValueError("Cannot specify dust_id when retrieving n_photons")

        # Check name
        available_components = self.get_available_components()
        if name not in available_components:
            raise Exception("name should be one of %s" % '/'.join(available_components))

        # If iteration is last one, find iteration number
        if iteration == -1:
            iteration = find_last_iteration(self.file)

        # Extract specific energy grid
        if name == 'temperature':
            array = np.array(self.file['iteration_%05i' % iteration]['specific_energy'])
            g_dust = self.file['Input/Dust']
            for i in range(array.shape[0]):
                dust = g_dust['dust_%03i' % (i + 1)]
                d = SphericalDust(dust)
                array[i, :, :, :] = d.mean_opacities._specific_energy2temperature(array[i, :, :, :])
        else:
            array = np.array(self.file['iteration_%05i' % iteration][name])

        # If required, extract grid for a specific dust type
        if name == 'n_photons':
            return array
        elif dust_id == 'all':
            return [array[i, :, :, :] for i in range(array.shape[0])]
        else:
            return array[dust_id, :, :, :]
