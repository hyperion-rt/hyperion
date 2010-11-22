import glob
import os

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

np.seterr(all='ignore')

from hyperion.dust import SimpleSphericalDust, prepare_emiss
from hyperion.util.integrate import interp1d_log10, integrate_loglog
from hyperion.util.constants import pi

# Read in Draine opacity file
d = SimpleSphericalDust('dustfiles_pah/draine_opac_new.dat')
d._extrapolate(1.e-3,1.e5)
d._sort()

# Create interpolating function for absorptive opacity
kappa_nu = interp1d_log10(d.nu, d.chi*(1.-d.albedo), bounds_error=False, fill_value=0.)

# Read in interstellar radiation field
isrf = np.loadtxt('dustfiles_pah/mmp83.dat', dtype=[('wav', float), ('jlam', float)])
wav = isrf['wav'][::-1] # wavelength in microns
jlam = isrf['jlam'][::-1] # 4*pi*J in ergs/cm^2/s/micron
nu = 3.e14 / wav

# Create interpolating function for mean intensity
jnu = interp1d_log10(nu, jlam * wav / nu, bounds_error=False, fill_value=0.)

# Find common frequency scale
nu_common = np.unique(np.hstack([d.nu, nu]))
nu_common.sort()

# Find the energy density of the interstellar radiation field weighted by the absorptive opacity
u_isrf = integrate_loglog(nu_common, jnu(nu_common) * kappa_nu(nu_common))

# Find the emissivity files
draine_emiss = glob.glob(os.path.join('dustfiles_pah/','emit_draine*.dat'))
U = [u_isrf * float(os.path.basename(filename).replace('emit_draine_', '').replace('.dat', '')) for filename in draine_emiss]

draine_emiss = np.array(draine_emiss)
U = np.array(U)
order = np.argsort(U)

prepare_emiss(draine_emiss[order], U[order], 'E', 'emit_draine.hdf5')

d.emissivities = 'emit_draine.hdf5'
d.write('hyperion/data/draine.hdf5')
d.plot('hyperion/data/draine.png')

