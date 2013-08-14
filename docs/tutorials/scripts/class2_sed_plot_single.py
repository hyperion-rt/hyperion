import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

# Open the model
m = ModelOutput('class2_sed.rtout')

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Extract the SED for the smallest inclination and largest aperture, and
# scale to 300pc. In Python, negative indices can be used for lists and
# arrays, and indicate the position from the end. So to get the SED in the
# largest aperture, we set aperture=-1.
sed = m.get_sed(inclination=0, aperture=-1, distance=300 * pc)

# Plot the SED. The loglog command is similar to plot, but automatically
# sets the x and y axes to be on a log scale.
ax.loglog(sed.wav, sed.val)

# Add some axis labels (we are using LaTeX here)
ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/s/cm$^2$]')

# Set view limits
ax.set_xlim(0.1, 2000.)
ax.set_ylim(2.e-16, 2.e-9)

# Write out the plot
fig.savefig('class2_sed_plot_single.png')
