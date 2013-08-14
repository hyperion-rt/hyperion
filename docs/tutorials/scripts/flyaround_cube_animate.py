import os

import numpy as np
import matplotlib.pyplot as plt

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

# Create output directory if it does not already exist
if not os.path.exists('frames'):
    os.mkdir('frames')

# Open model
m = ModelOutput('flyaround_cube.rtout')

# Read image from model
image = m.get_image(distance=300 * pc, units='MJy/sr')

# image.val is now an array with four dimensions (n_view, n_y, n_x, n_wav)

for iview in range(image.val.shape[0]):

    # Open figure and create axes
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(1, 1, 1)

    # This is the command to show the image. The parameters vmin and vmax are
    # the min and max levels for the grayscale (remove for default values).
    # The colormap is set here to be a heat map. Other possible heat maps
    # include plt.cm.gray (grayscale), plt.cm.gist_yarg (inverted grayscale),
    # plt.cm.jet (default, colorful). The np.sqrt() is used to plot the
    # images on a sqrt stretch.
    ax.imshow(np.sqrt(image.val[iview, :, :, 0]), vmin=0, vmax=np.sqrt(2000.),
              cmap=plt.cm.gist_heat, origin='lower')

    # Save figure. The facecolor='black' and edgecolor='black' are for
    # esthetics, and hide the axes
    fig.savefig('frames/frame_%05i.png' % iview,
                facecolor='black', edgecolor='black')

    # Close figure
    plt.close(fig)
