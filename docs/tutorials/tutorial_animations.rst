.. _animations:

=================
Making animations
=================

The tutorial here assumes that you are using the :ref:`tutorial-model`.

The following script describes how to generate PNG frames for an animation::

    import os

    import numpy as np

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc

    # Create output directory if it does not already exist
    if not os.path.exists('frames'):
        os.mkdir('frames')

    # Open model
    m = ModelOutput('tutorial_model.rtout')

    # Read image from model
    wav, nufnu = m.get_image(group=2, distance=300 * pc)

    # nufnu is now an array with four dimensions (n_view, n_wav, n_y, n_x)

    # Fix the wavelength to the first one and cycle through viewing angles
    iwav = 0
    print "Wavelength is %g microns" % wav[iwav]

    for iview in range(nufnu.shape[0]):

        # Open figure and create axes
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # This is the command to show the image. The parameters vmin and vmax are
        # the min and max levels for the grayscale (remove for default values).
        # The colormap is set here to be a heat map. Other possible heat maps
        # include plt.cm.gray (grayscale), plt.cm.gist_yarg (inverted grayscale),
        # plt.cm.jet (default, colorful). The np.sqrt() is used to plot the
        # images on a sqrt stretch.
        ax.imshow(np.sqrt(nufnu[iview, :, :, iwav]), vmin=0, vmax=np.sqrt(1.e-11), \
                  cmap=plt.cm.gist_heat, origin='lower')

        # Save figure. The facecolor='black' and edgecolor='black' are for
        # esthetics, and hide the axes
        fig.savefig('frames/frame_%05i.png' % iview, \
                    facecolor='black', edgecolor='black')

        # Close figure
        plt.close(fig)

The frames can then be combined into a GIF animation using ImageMagick::

    $ convert -delay 10 -adjoin frames/*.png movie.gif

The delay value is the delay between frames in 1/100ths of a second. The result is the following:

.. image:: images/movie.gif
   :scale: 75 %
   :alt: Fly-around movie of simple model
   :align: center