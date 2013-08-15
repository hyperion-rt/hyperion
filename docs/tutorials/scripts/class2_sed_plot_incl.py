import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

m = ModelOutput('class2_sed.rtout')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Extract all SEDs
sed = m.get_sed(inclination='all', aperture=-1, distance=300 * pc)

# Plot SED for each inclination
for i in range(sed.val.shape[0]):
    ax.loglog(sed.wav, sed.val[i, :], color='black')

ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/s/cm$^2$]')
ax.set_xlim(0.1, 2000.)
ax.set_ylim(2.e-16, 2.e-9)
fig.savefig('class2_sed_plot_incl.png')
