import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from hyperion.util.constants import pc

mo = ModelOutput('pure_scattering.rtout')

image_fnu = mo.get_image(inclination=0, units='MJy/sr', distance=300. * pc)
image_pol = mo.get_image(inclination=0, stokes='linpol')

fig = plt.figure(figsize=(8, 8))

# Make total intensity sub-plot

ax = fig.add_axes([0.1, 0.3, 0.4, 0.4])
ax.imshow(image_fnu.val[:, :, 0], extent=[-13, 13, -13, 13],
          interpolation='none', cmap=plt.cm.gist_heat,
          origin='lower', vmin=0., vmax=4e9)
ax.set_xlim(-13., 13.)
ax.set_ylim(-13., 13.)
ax.set_xlabel("x (solar radii)")
ax.set_ylabel("y (solar radii)")
ax.set_title("Surface brightness")

# Make linear polarization sub-plot

ax = fig.add_axes([0.51, 0.3, 0.4, 0.4])
im = ax.imshow(image_pol.val[:, :, 0] * 100., extent=[-13, 13, -13, 13],
               interpolation='none', cmap=plt.cm.gist_heat,
               origin='lower', vmin=0., vmax=100.)
ax.set_xlim(-13., 13.)
ax.set_ylim(-13., 13.)
ax.set_xlabel("x (solar radii)")
ax.set_title("Linear Polarization")
ax.set_yticklabels('')

axcb = fig.add_axes([0.92, 0.3, 0.02, 0.4])
plt.colorbar(im, label="%", cax=axcb)
fig.savefig('pure_scattering_inner_disk.png', bbox_inches='tight')
