from __future__ import print_function, division

import os
import warnings

import h5py
import numpy as np

from ..util.constants import c, pi
from ..util.functions import FreezableClass
from ..dust import SphericalDust
from astropy import log as logger
from ..util.decorator import decorator
from ..grid import CartesianGrid, SphericalPolarGrid, CylindricalPolarGrid, OctreeGrid, AMRGrid

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
    '''
    A class that can be used to access data in the output file from
    radiative transfer models.

    Parameters
    ----------
    name : str
        The name of the model output file (including extension)
    '''

    def __init__(self, filename):

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
        Retrieve SEDs for a specific image group and Stokes component.

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
            The index of the aperture to plot (zero-based). Use 'all' to
            return all apertures, and -1 to show the largest aperture.

        inclination : int, optional
            The index of the viewing angle to plot (zero-based). Use 'all'
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
            specific source and dust component can be specified (where 0 is
            the first source or dust type). If 'all' is specified, then all
            components are returned individually. If neither of these are
            not specified, then the total component requested for all
            sources or dust types is returned. For sources, it is also possible
            to specify a source name as a string, if the source name was set
            during the set-up.

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

        uncertainty : numpy.ndarray
            The uncertainties on the flux or degree of polarization. This
            has the same dimensions as the flux or degree of polarization
            array. This is only returned if uncertainties were requested.
        '''

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

                    if isinstance(source_id, basestring) and source_id != 'all':
                        try:
                            source_id = self.get_available_sources().index(source_id)
                        except ValueError:
                            raise ValueError("No source named {0}".format(source_id))

                    ns = g['seds'].attrs['n_sources']
                    nd = g['seds'].attrs['n_dust']

                    io = ((io - (io + 1) % 2 + 1) * ns + (io - io % 2) * nd) / 2

                    if component.startswith('source'):
                        if source_id is None or source_id == 'all':
                            io = (io, io + ns)
                        else:
                            if source_id < 0 or source_id >= ns:
                                raise ValueError("source_id should be between 0 and %i" % (ns - 1))
                            io = io + source_id
                    else:
                        if dust_id is None or dust_id == 'all':
                            io = (io, io + nd)
                        else:
                            if dust_id < 0 or dust_id >= nd:
                                raise ValueError("dust_id should be between 0 and %i" % (nd - 1))
                            io = io + dust_id

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

        # Set default units
        if units is None:
            if distance is None and not inside_observer:
                units = 'ergs/s'
            else:
                units = 'ergs/cm^2/s'

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
                raise ValueError("Unknown units: %s" % units)

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
            if not isinstance(aperture, int):
                raise TypeError('aperture should be an integer (it should'
                                ' be the index of the aperture, not the '
                                'value itself)')
            flux = flux[:, :, :, aperture]
            if uncertainties:
                unc = unc[:, :, :, aperture]

        if inclination == 'all':
            pass
        else:
            if not isinstance(inclination, int):
                raise TypeError('inclination should be an integer (it should'
                                ' be the index of the inclination, not the '
                                'value itself)')
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

        from .sed import SED
        sed = SED(nu=nu, val=flux, unc=unc if uncertainties else None, units=units)

        # Add aperture information
        sed.ap_min = g['seds'].attrs['apmin']
        sed.ap_max = g['seds'].attrs['apmax']

        # Add distance
        sed.distance = distance

        # Save whether the SED was from an inside observer
        sed.inside_observer = inside_observer

        return sed

    @on_the_fly_hdf5
    def get_image(self, stokes='I', group=0, technique='peeled',
                  distance=None, component='total', inclination='all',
                  uncertainties=False, units=None,
                  source_id=None, dust_id=None):
        '''
        Retrieve images for a specific image group and Stokes component.

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
            The index of the viewing angle to plot (zero-based). Use 'all'
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
            specific source and dust component can be specified (where 0 is
            the first source or dust type). If 'all' is specified, then all
            components are returned individually. If neither of these are
            not specified, then the total component requested for all
            sources or dust types is returned. For sources, it is also possible
            to specify a source name as a string, if the source name was set
            during the set-up.

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

        uncertainty : numpy.ndarray
            The uncertainties on the flux or degree of polarization. This
            has the same dimensions as the flux or degree of polarization
            array. This is only returned if uncertainties were requested.
        '''

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

            track_origin = g['images'].attrs['track_origin'].decode('utf-8')

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

                    if isinstance(source_id, basestring) and source_id != 'all':
                        try:
                            source_id = self.get_available_sources().index(source_id)
                        except ValueError:
                            raise ValueError("No source named {0}".format(source_id))

                    ns = g['images'].attrs['n_sources']
                    nd = g['images'].attrs['n_dust']

                    io = ((io - (io + 1) % 2 + 1) * ns + (io - io % 2) * nd) / 2

                    if component.startswith('source'):
                        if source_id is None or source_id == 'all':
                            io = (io, io + ns)
                        else:
                            if source_id < 0 or source_id >= ns:
                                raise ValueError("source_id should be between 0 and %i" % (ns - 1))
                            io = io + source_id
                    else:
                        if dust_id is None or dust_id == 'all':
                            io = (io, io + nd)
                        else:
                            if dust_id < 0 or dust_id >= nd:
                                raise ValueError("dust_id should be between 0 and %i" % (nd - 1))
                            io = io + dust_id

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

        # Set default units
        if units is None:
            if distance is None and not inside_observer:
                units = 'ergs/s'
            else:
                units = 'ergs/cm^2/s'

        # Find pixel dimensions of image
        ny, nx = flux.shape[-3:-1]

        # Find physical and angular extent, and pixel area in steradians
        if inside_observer:

            # Physical extent cannot be set
            x_min = x_max = y_min = y_max = None

            # Find extent of the image
            lon_min = g['images'].attrs['xmin']
            lon_max = g['images'].attrs['xmax']
            lat_min = g['images'].attrs['ymin']
            lat_max = g['images'].attrs['ymax']

            # Need to construct a mesh since every pixel might have a
            # different size
            lon = np.linspace(np.radians(lon_min), np.radians(lon_max), nx + 1)
            lat = np.cos(np.linspace(np.radians(90. - lat_min), np.radians(90. - lat_max), ny + 1))
            dlon = lon[1:] - lon[:-1]
            dlat = lat[:-1] - lat[1:]
            DLON, DLAT = np.meshgrid(dlon, dlat)

            # Find pixel area in steradians
            pix_area_sr = DLON * DLAT

        else:

            # Find physical extent of the image
            x_min = g['images'].attrs['xmin']
            x_max = g['images'].attrs['xmax']
            y_min = g['images'].attrs['ymin']
            y_max = g['images'].attrs['ymax']

            if distance is not None:

                # Find angular extent
                lon_min_rad = np.arctan(x_min / distance)
                lon_max_rad = np.arctan(x_max / distance)
                lat_min_rad = np.arctan(y_min / distance)
                lat_max_rad = np.arctan(y_max / distance)

                # Find pixel size in arcseconds
                pix_dx = abs(lon_max_rad - lon_min_rad) / float(nx)
                pix_dy = abs(lat_max_rad - lat_min_rad) / float(ny)

                # Find pixel area in steradians
                pix_area_sr = pix_dx * pix_dy

                # Find angular extent in degrees
                lon_min = np.degrees(lon_min_rad)
                lon_max = np.degrees(lon_max_rad)
                lat_min = np.degrees(lat_min_rad)
                lat_max = np.degrees(lat_max_rad)

            else:

                # Angular extent cannot be set
                lon_min = lon_max = lat_min = lat_max = None

                # Pixel area in steradians cannot be set
                pix_area_sr = None

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
            elif units == 'MJy/sr':
                if inside_observer:
                    scale = 1.e17 / nu[np.newaxis, np.newaxis, :] / pix_area_sr[:, :, np.newaxis]
                else:
                    scale = 1.e17 / nu / pix_area_sr
            else:
                raise ValueError("Unknown units: %s" % units)

            # Scale by distance
            if distance:
                scale *= 1. / (4. * pi * distance ** 2)

        else:

            if units != 'ergs/s':
                raise ValueError("Since distance= is not specified, units should be set to ergs/s")

            scale = np.repeat(1., len(nu))

        # If in 32-bit mode, need to convert to 64-bit because of
        # scaling/polarization to be safe
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
            if not isinstance(inclination, int):
                raise TypeError('inclination should be an integer (it should'
                                ' be the index of the inclination, not the '
                                'value itself)')
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

        from .image import Image
        image = Image(nu=nu, val=flux, unc=unc if uncertainties else None, units=units)

        # Add physical extent
        image.x_min = x_min
        image.x_max = x_max
        image.y_min = y_min
        image.y_max = y_max

        # Add angular extent
        image.lon_min = lon_min
        image.lon_max = lon_max
        image.lat_min = lat_min
        image.lat_max = lat_max

        # Add pixel area in steradians
        image.pix_area_sr = pix_area_sr

        # Add distance
        image.distance = distance

        # Save whether the image was from an inside observer
        image.inside_observer = inside_observer

        return image

    @on_the_fly_hdf5
    def get_available_sources(self):
        """
        Find out what sources are available in the output image (useful if detailed tracking was used)
        """
        if self.file['Input'].file != self.file.file:
            # Workaround for h5py bug - can't access link directly,
            # need to use file attribute
            g_sources = self.file['Input'].file[self.file['Input'].name]['Sources']
        else:
            g_sources = self.file['Input/Sources']

        return sorted([s.attrs['name'].decode('utf-8') for s in g_sources.values()])

    @on_the_fly_hdf5
    def get_available_components(self, iteration=-1):
        '''
        Find out what physical components are available in the output file

        Parameters
        ----------

        iteration : integer, optional
            The iteration to retrieve the grid for. The default is to return the components for the last iteration
        '''

        from .helpers import find_last_iteration

        # If iteration is last one, find iteration number
        if iteration == -1:
            iteration = find_last_iteration(self.file)

        # Return components
        components = list(self.file['iteration_%05i' % iteration].keys())
        if 'specific_energy' in components:
            components.append('temperature')
        return components

    @on_the_fly_hdf5
    def get_quantities(self, iteration=-1):
        '''
        Retrieve one of the physical grids for the model

        Parameters
        ----------
        iteration : integer, optional
            The iteration to retrieve the quantities for. The default is to return the grid for the last iteration.

        Returns
        -------

        grid : Grid instance
            An object containing information about the geometry and quantities
        '''

        from .helpers import find_last_iteration


        if self.file['Input'].file != self.file.file:
            # Workaround for h5py bug - can't access link directly,
            # need to use file attribute
            g_grid = self.file['Input'].file[self.file['Input'].name]['Grid']
        else:
            g_grid = self.file['Input/Grid']

        # Find coordinate grid type
        coord_type = g_grid['Geometry'].attrs['grid_type'].decode('utf-8')

        if coord_type == 'car':
            g = CartesianGrid()
        elif coord_type == 'cyl_pol':
            g = CylindricalPolarGrid()
        elif coord_type == 'sph_pol':
            g = SphericalPolarGrid()
        elif coord_type == 'amr':
            g = AMRGrid()
        elif coord_type == 'oct':
            g = OctreeGrid()

        # Read in geometry and input quantities
        g.read_geometry(g_grid['Geometry'])

        # If iteration is last one, find iteration number
        if iteration == -1:
            iteration = find_last_iteration(self.file)
        else:
            iteration = iteration + 1  # Python value is zero based

        # Read in quantities from the requested iteration (if available)
        if iteration > 0:
            g.read_quantities(self.file['iteration_%05i' % iteration])

        if not 'density' in g:
            logger.info("No density present in output, reading initial density")
            g.read_quantities(g_grid['Quantities'], quantities=['density'])

        # Compute the temperature as a derived quantity
        if 'specific_energy' in g:

            # Open the dust group
            n_dust = g.n_dust
            g_dust = self.file['Input/Dust']

            # Compile a list of specific energy to temperature functions
            convert_func = []
            for i in range(n_dust):

                # Read in dust type
                dust = g_dust['dust_%03i' % (i + 1)]
                d = SphericalDust(dust)

                # Add to conversion function list
                convert_func.append(d.mean_opacities._specific_energy2temperature)

            # Define function to convert different specific energy to
            # temperature for different dust types
            def specific_energy2temperature(quantities):
                quantities['temperature'] = []
                for i in range(n_dust):
                    quantities['temperature'].append(convert_func[i](quantities['specific_energy'][i]))

            # Get the grid to add the quantity
            g.add_derived_quantity('temperature', specific_energy2temperature)

        return g

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

        from .helpers import find_last_iteration

        warnings.warn("get_physical_grid is deprecated, use get_quantities instead", DeprecationWarning)

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
            g_dust = g_dust.file[g_dust.name]  # workaround for h5py < 2.1.0
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
