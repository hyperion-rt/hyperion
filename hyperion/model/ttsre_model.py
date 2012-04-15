from __future__ import print_function, division

import os
import shutil

import numpy as np
from atpy import Table

import hyperion
from ..util.parfile import parse
from .. import atmos
from ..dust import SimpleSphericalDust
from . import AnalyticalYSOModel
from ..util.functions import filename2hdf5
from ..util.constants import msun, rsun, au, year, k, m_h, G, pi, sigma, c
from ..util.logger import logger


class TtsreModel(AnalyticalYSOModel):

    def load_from_mctherm(self, filename):

        tsub = 1600.

        # Read in mctherm.par file
        par = parse(filename)

        if not 'nrg' in par:
            par['nrg'] = 400

        if not 'ntg' in par:
            par['ntg'] = 300

        if not 'npg' in par:
            par['npg'] = 2

        # Remove one from number of grid cells (actually number of walls)
        self.nr = par['nrg'] - 1
        self.nt = par['ntg'] - 1
        self.np = par['npg'] - 1

        if not 'fsg(1)' in par:
            par['fsg(1)'] = 0.

        if not 'fsg(2)' in par:
            par['fsg(2)'] = 0.

        if not 'fsg(3)' in par:
            par['fsg(3)'] = 0.

        if not 'fsg(4)' in par:
            par['fsg(4)'] = 0.

        self.disk1_pah_frac = par['fsg(1)']
        self.disk2_pah_frac = par['fsg(2)']
        self.envelope_pah_frac = par['fsg(3)']
        self.cavity_pah_frac = par['fsg(4)']

        self.disk1.cylindrical_inner_rim = False
        self.disk1.cylindrical_outer_rim = False
        self.disk2.cylindrical_inner_rim = False
        self.disk2.cylindrical_outer_rim = False

        # STELLAR PARAMETERS

        mstar = par['massc'] * msun
        tstar = par['tstar']
        rstar = par['rstar'] * rsun
        lstar = 4. * pi * rstar ** 2 * sigma * tstar ** 4

        self.set_stellar_radius(rstar)
        self.star.luminosity = lstar
        self.star.limb = par['climb']

        # STELLAR SPECTRUM

        self.planck = par['cplanckst']

        if not 'pardir' in par:
            par['pardir'] = './'

        if not self.planck:

            atmos_file_orig = os.path.abspath(os.path.dirname(filename)) + "/" + par['pardir'] + "/" + par['atname']
            atmos_file = filename2hdf5(atmos_file_orig)

            # Check if file already exists in database
            try:
                path = hyperion.get_HDF5_datafile(atmos_file_orig)
                logger.info("Exact HDF5 match found for %s" % os.path.basename(atmos_file_orig))
                atmos_file = Table(path)
            except:
                logger.info("No existing HDF5 file for %s" % os.path.basename(atmos_file_orig))
                logger.info(" -> computing from scratch (may take a few minutes)")
                atmos_file = atmos.prepare_atmos(atmos_file_orig)

            self.set_stellar_spectrum(atmos_file)

        else:

            self.set_stellar_temperature(tstar)

        # DUST FILES

        if not 'dustname(5)' in par:
            par['dustname(5)'] = par['dustname(1)']

        if not 'dustname(6)' in par:
            par['dustname(6)'] = par['dustname(2)']

        if not 'dustname(7)' in par:
            par['dustname(7)'] = par['dustname(3)']

        if not 'dustname(8)' in par:
            par['dustname(8)'] = par['dustname(4)']

        for i, parameter in enumerate(['disk1_dust', 'disk2_dust',
                                       'envelope_dust', 'cavity_dust',
                                       'disk1_pah', 'disk2_pah',
                                       'envelope_pah', 'cavity_pah']):

            # Read in dust file
            dust_file_orig = os.path.abspath(os.path.dirname(filename)) + "/" + par['pardir'] + "/" + par['dustname(%i)' % (i + 1)]
            dust_file = filename2hdf5(dust_file_orig)

            # Check if file already exists in database
            try:
                path = hyperion.get_HDF5_datafile(dust_file_orig)
                logger.info("Exact HDF5 match found for %s" % os.path.basename(dust_file_orig))
                shutil.copy(path, dust_file)
            except:
                logger.info("No existing HDF5 file for %s" % os.path.basename(dust_file_orig))
                logger.info(" -> computing from scratch (may take a few minutes)")
                d = SimpleSphericalDust(dust_file_orig)
                d._extrapolate(1.e-3, 1.e5)
                d.write(dust_file)

            self.__dict__[parameter] = dust_file

        # DISK PROPERTIES

        if not 'fmassd1' in par:
            par['fmassd1'] = 0.

        if not 'a(1)' in par:
            par['a(1)'] = par['a']

        if not 'b(1)' in par:
            par['b(1)'] = par['b']

        if not 'a(2)' in par:
            par['a(2)'] = par['a']

        if not 'b(2)' in par:
            par['b(2)'] = par['b']

        if not 'zscale(1)' in par:
            par['zscale(1)'] = par['zscale']

        if not 'zscale(2)' in par:
            par['zscale(2)'] = par['zscale']

        self.disk1.mstar = mstar
        self.disk2.mstar = mstar
        self.disk1.mass = par['massd'] * msun * par['fmassd1']
        self.disk2.mass = par['massd'] * msun * (1. - par['fmassd1'])
        self.disk1.rmax = par['rmaxd'] * au
        self.disk2.rmax = par['rmaxd'] * au
        self.disk1.alpha = par['a(1)']
        self.disk2.alpha = par['a(2)']
        self.disk1.beta = par['b(1)']
        self.disk2.beta = par['b(2)']

        rtrunc = par['rtrunc'] * rstar

        # ENVELOPE PROPERTIES

        self.envelope.mstar = mstar
        self.envelope.mdot = par['rate'] * msun / year
        self.envelope.rmax = par['rmax'] * au
        self.envelope.rc = par['rc'] * au
        self.envelope.rchole = par['rchole'] * au
        self.envelope.rhoamb = par['rhoamb']

        # BIPILAR CAVITY

        self.cavity.shape = par['cshape'].lower()
        self.cavity.theta = par['thet1']
        self.cavity.exp = par['ex1']
        self.cavity.rmax = self.envelope.rmax
        self.cavity.exf = par['exf']
        self.cavity.rhoconst = par['rhoconst1']
        self.cavity.rhoamb = par['rhoamb']

        # Compute dust sublimation radius
        rsub = rstar * (tstar / tsub) ** 2.1

        # DISK INNER RADII

        if par['crminsub']:
            self.disk1.rmin = max(rtrunc, par['rmind'] * rsub)
            self.disk2.rmin = max(rtrunc, par['rmind'] * rsub)
            self.cavity.rmin = max(rtrunc, par['rmine'] * rsub)
            self.envelope.rmin = max(rtrunc, par['rmine'] * rsub)
            # IMPLEMENT: rmine2
            # IMPLEMENT: rmin_sg
        else:
            self.disk1.rmin = max(rtrunc, par['rmind'] * rstar)
            self.disk2.rmin = max(rtrunc, par['rmind'] * rstar)
            self.cavity.rmin = max(rtrunc, par['rmine'] * rstar)
            self.envelope.rmin = max(rtrunc, par['rmine'] * rstar)
            # IMPLEMENT: rmine2
            # IMPLEMENT: rmin_sg

        # DISK SCALEHEIGHTS

        if par['czmin'] is True:
            par['czmin'] = 'rsub'

        if par['czmin'] is False:
            par['czmin'] = 'rstar'

        if par['czmin'].lower() == 'rstar':

            # Scale by R_star
            self.disk1.h_star = par['zscale(1)'] * rstar
            self.disk2.h_star = par['zscale(2)'] * rstar

        elif par['czmin'].lower() == 'rsub':

            # Scale by hydrostatic scaleheight
            c_s = np.sqrt(k * tsub / 2.3 / m_h)
            z_hseq_rsub = rsub * c_s * np.sqrt(rsub / G / mstar)
            self.disk1.h_star = par['zscale(1)'] * z_hseq_rsub * (rstar / rsub) ** par['b(1)']
            self.disk2.h_star = par['zscale(2)'] * z_hseq_rsub * (rstar / rsub) ** par['b(2)']

        elif par['czmin'].lower() == 'r100':

            # Scale by AU
            self.disk1.h_star = par['zscale(1)'] * au * (rstar / (100. * au)) ** par['b(1)']
            self.disk2.h_star = par['zscale(2)'] * au * (rstar / (100. * au)) ** par['b(2)']

        else:

            raise Exception("Unknown CZMIN: %s" % par['czmin'])

        # DISK ACCRETION RATE

        self.accretion = par['cdiskacc']

        if self.accretion:

            if not par['calpha']:
                mdotdisk = par['alpha'] * msun / year
            else:
                mdotdisk1 = np.sqrt(18 * pi ** 3) * par['alpha'] * np.sqrt(G * mstar / rstar) * self.disk1.rho_0() * self.disk1.h_star ** 3 / rstar
                mdotdisk2 = np.sqrt(18 * pi ** 3) * par['alpha'] * np.sqrt(G * mstar / rstar) * self.disk2.rho_0() * self.disk2.h_star ** 3 / rstar
                mdotdisk = mdotdisk1 + mdotdisk2

            self.disk1.mdot = mdotdisk * par['fmassd1']
            self.disk2.mdot = mdotdisk * (1. - par['fmassd1'])

            # STELLAR ACCRETION LUMINOSITY

            lshock = G * mstar * mdotdisk * (1 / rstar - 1 / rtrunc)  # W

            if not 'cspot' in par:
                par['cspot'] = False

            if par['cspot']:

                if par['nspot'] == 1:
                    spotsize = np.degrees(np.arccos(1. - 2. * par['fspot']))
                elif par['nspot'] == 2:
                    spotsize = np.degrees(np.arccos(1. - par['fspot']))
                else:
                    raise Exception("Invalid NSPOT: %s" % par['nspot'])

                fluxratio = 0.5 * lshock / lstar / par['fspot']

            else:

                fluxratio = 0.5 * lshock / lstar

            tshock = tstar * (1 + fluxratio) ** 0.25  # K

            # Create x-ray spectrum file
            wav = np.linspace(0.015, 0.060, 100)
            xray = Table()
            xray.add_column('nu', c * 1.e4 / wav)
            xray.add_column('wav', wav)
            xray.add_column('fnu', 1.)
            xray.sort('nu')

            if par['cspot']:

                if par['nspot'] == 1:
                    self.add_spot(lshock / 2., longitude=0., latitude=par['spotlat'], size=spotsize, spectrum=xray)
                    self.add_spot(lshock / 2., longitude=0., latitude=par['spotlat'], size=spotsize, temperature=tshock)
                else:
                    self.add_spot(lshock / 4., longitude=0., latitude=par['spotlat'], size=spotsize, spectrum=xray)
                    self.add_spot(lshock / 4., longitude=0, latitude=par['spotlat'], size=spotsize, temperature=tshock)
                    self.add_spot(lshock / 4., longitude=180., latitude=90. - par['spotlat'], size=spotsize, spectrum=xray)
                    self.add_spot(lshock / 4., longitude=180., latitude=90. - par['spotlat'], size=spotsize, temperature=tshock)

            else:

                self.star_uv.set_luminosity(lshock / 2.)
                self.star_uv.set_temperature(tshock)
                self.star_uv.limb = par['climb']

                self.star_xray.set_luminosity(lshock / 2.)
                self.star_xray.set_spectrum(xray)
                self.star_xray.limb = par['climb']

        # IMAGES/SEDS

        if not 'imcube' in par:
            par['imcube'] = False

        if not 'nfreq' in par:
            par['nfreq'] = 250

        if not 'phie' in par:
            par['phie'] = 30.

        if par['imcube']:
            imsize = 200
        else:
            imsize = 1

        nwav = par['nfreq']
        wav_min = 0.01
        wav_max = 5000.

        par['rmaxi'] *= au
        par['apmin'] *= au
        par['apmax'] *= au

        image = self.add_binned_images()
        image.set_viewing_bins(10, 1)
        image.set_image_size(imsize, imsize)
        image.set_image_limits(-par['rmaxi'], par['rmaxi'], -par['rmaxi'], par['rmaxi'])
        image.set_aperture_range(par['nap'], par['apmin'], par['apmax'])
        image.set_wavelength_range(nwav, wav_min, wav_max)

        if par['cpeel']:

            if type(par['thete']) == str:
                theta = np.array(par['thete'].split(), dtype=float)
                phi = np.array(par['phie'].split(), dtype=float)
            else:
                theta = np.array([par['thete']], dtype=float)
                phi = np.array([par['phie']], dtype=float)

            image = self.add_peeled_images()
            image.set_viewing_angles(theta, phi)
            image.set_image_size(imsize, imsize)
            image.set_image_limits(-par['rmaxi'], par['rmaxi'], -par['rmaxi'], par['rmaxi'])
            image.set_aperture_range(par['nap'], par['apmin'], par['apmax'])
            image.set_wavelength_range(nwav, wav_min, wav_max)

        # OVERALL CONFIGURATION

        if not 'diffusion' in par:
            par['diffusion'] = True

        if not 'nimin' in par:
            par['nimin'] = 5

        if not 'npmin' in par:
            par['npmin'] = par['np']

        self.set_pda(par['diffusion'])
        self.set_mrw(par['diffusion'])
        self.set_n_lucy_iterations(max(par['nimin'], 5))
        self.set_n_lucy_photons(par['npmin'])
        self.set_n_last_photons(par['np'])
        self.set_n_ray_photons(par['np'])
        self.set_n_stats(par['iwrite'])

        self.set_spherical_polar_grid_auto(self.nr, self.nt, self.np)
