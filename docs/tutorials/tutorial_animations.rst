.. _animations:

=================
Making animations
=================

In this tutorial, we will find out how to make animations from Hyperion output. To start with we will use a model very similar to :doc:`tutorial_images` but this time we only compute the image in one wavelength bin, but for a number of different viewing angles.

The tutorial here assumes that you are using the :ref:`tutorial-model`.

The following script describes how to generate PNG frames for an animation:

.. literalinclude:: scripts/flyaround_cube_animate.py
   :language: python

The frames can then be combined into a GIF animation using ImageMagick::

    $ convert -delay 10 -adjoin frames/*.png movie.gif

The delay value is the delay between frames in 1/100ths of a second. The result is the following:

.. image:: scripts/movie.gif
   :alt: Fly-around movie of simple model
   :align: center