import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

mo = ModelOutput('class1_example.rtout')
sed = mo.get_sed(aperture=-1, distance=140. * pc)

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(1, 1, 1)
ax.loglog(sed.wav, sed.val.transpose(), color='black')
ax.set_xlim(0.03, 2000.)
ax.set_ylim(2.e-15, 1e-8)
ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/cm$^2/s$]')
fig.savefig('class1_example_sed.png', bbox_inches='tight')
