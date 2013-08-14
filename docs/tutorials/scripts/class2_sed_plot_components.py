import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

m = ModelOutput('class2_sed.rtout')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Total SED
sed = m.get_sed(inclination=0, aperture=-1, distance=300 * pc)
ax.loglog(sed.wav, sed.val, color='black', lw=3, alpha=0.5)

# Direct stellar photons
sed = m.get_sed(inclination=0, aperture=-1, distance=300 * pc,
                       component='source_emit')
ax.loglog(sed.wav, sed.val, color='blue')

# Scattered stellar photons
sed = m.get_sed(inclination=0, aperture=-1, distance=300 * pc,
                       component='source_scat')
ax.loglog(sed.wav, sed.val, color='teal')

# Direct dust photons
sed = m.get_sed(inclination=0, aperture=-1, distance=300 * pc,
                       component='dust_emit')
ax.loglog(sed.wav, sed.val, color='red')

# Scattered dust photons
sed = m.get_sed(inclination=0, aperture=-1, distance=300 * pc,
                       component='dust_scat')
ax.loglog(sed.wav, sed.val, color='orange')


ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/s/cm$^2$]')
ax.set_xlim(0.1, 2000.)
ax.set_ylim(2.e-16, 2.e-9)
fig.savefig('class2_sed_plot_components.png')
