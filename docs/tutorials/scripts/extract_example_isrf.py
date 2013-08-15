import numpy as np
import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.integrate import integrate_loglog

# Use LaTeX for plots
plt.rc('text', usetex=True)

# Open the output file
m = ModelOutput('example_isrf.rtout')

# Get an all-sky flux map
image = m.get_image(units='ergs/cm^2/s/Hz', inclination=0)

# Compute the frequency-integrated flux
fint = np.zeros(image.val.shape[:-1])
for (j, i) in np.ndindex(fint.shape):
    fint[j, i] = integrate_loglog(image.nu, image.val[j, i, :])

# Find the area of each pixel
l = np.radians(np.linspace(180., -180., fint.shape[1] + 1))
b = np.radians(np.linspace(-90., 90., fint.shape[0] + 1))
dl = l[1:] - l[:-1]
db = np.sin(b[1:]) - np.sin(b[:-1])
DL, DB = np.meshgrid(dl, db)
area = np.abs(DL * DB)

# Compute the intensity
intensity = fint / area

# Intitialize plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='aitoff')

# Show intensity
image = ax.pcolormesh(l, b, intensity, cmap=plt.cm.gist_heat, vmin=0.0, vmax=0.0025)

# Add mean intensity
four_pi_jnu = round(np.sum(intensity * area), 4)
fig.text(0.40, 0.15, r"$4\pi J = %6.4f$ "
                     r"${\rm erg\,cm^{-2}\,s^{-1}}$" % four_pi_jnu, size=14)

# Add a colorbar
cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
fig.colorbar(image, cax=cax)

# Add title and improve esthetics
ax.set_title(r"Integrated intensity in a given pixel "
             r"(${\rm erg\,cm^{-2}\,s^{-1}\,ster^{-1}}$)", size=12, y=1.1)
ax.grid()
ax.tick_params(axis='both', which='major', labelsize=10)
cax.tick_params(axis='both', which='major', labelsize=10)

# Save the plot
fig.savefig('isrf_intensity.png', bbox_inches='tight')
