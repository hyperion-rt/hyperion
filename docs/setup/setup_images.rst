Setting up images and SEDs
==========================

There are two main kinds of images/SEDs that can be produced for each model:
images/SEDs computed by binning the photons as they escape from the density
grid, and images/SEDs computed by peeling off photon packets at each
interaction into well defined directions. The latter provide more accurate
SEDs and much better signal-to-noise, and are likely to be more commonly used
than the former.

The code currently allows at most one set of binned images, and any number
of sets of peeled images. A set is defined by a wavelength range, image
resolution and extent, and any number of viewing angles.

Common parameters
-----------------

The wavelength range (in microns) for the images/SEDs can be specified using::

    image.set_wavelength_range(n_wav, wav_min, wav_max)

The image size in pixels and the extent of the images can be specified using::

    image.set_image_size(n_x, n_y)
    image.set_image_limits(xmin, xmax, ymin, ymax)

where the image limits should be given in cm. The apertures for the SEDs can
be specified using::

    image.set_aperture_range(n_ap, ap_min, ap_max)

where the apertures should be given in cm. The default is to have one
aperture with infinite size, i.e. measuring all the flux.

Binned images
-------------

To add a set of binned images/SEDs to the model, use::

    image = m.add_binned_images()

The number of bins in the theta and phi direction can be specified using::

    image.set_viewing_bins(10, 10)

Peeled images
-------------

To add a set of peeled images/SEDs to the model, use::

    image = m.add_peeled_images()

The viewing angles should be specified as lists or arrays of theta and phi
values, in degrees. For example, the following produces images from pole-on
to edge-on at constant phi using 20 viewing angles::

    # Set number of viewing angles
    n_view = 20

    # Generate the viewing angles
    theta = np.linspace(0., 90., n_view)
    phi = np.repeat(45., n_view)

    # Set the viewing angles
    image.set_viewing_angles(theta, phi)

A few more advanced parameters are also available, and these are described in :doc:`../advanced/peeloff`.

.. note:: For peeled images, the number of viewing angles directly impacts the
          performance of the code - once the specific energy/temperature has
          been computed, the code will then run approximately in a time
          proportional to the number of viewing angles.

Uncertainties
-------------

Uncertainties can be computed for SEDs/images (doubling the memory/disk space required)::

    image.set_uncertainties(True)

File output
-----------

Finally, to save space, images can be written out as 32-bit floats instead of 64-bit floats. To write them out as 32-bit floats, use::

    image.set_output_bytes(4)

and to write them out as 64-bit floats, use::

    image.set_output_bytes(8)

Tracking photon origin
----------------------

SEDs/images can also be split into emitted/thermal or scattered components
from sources or dust (4 combinations). To activate this, use::

    image.set_track_origin('basic')

It is also possible to split the SED into individual sources and dust types::

    image.set_track_origin('detailed')

For example, if five sources and two dust types are present, there will be 14
components in total: five for photons emitted from sources, two for photons
emitted from dust, five for photons emitted from sources and subsequently
scattered, and two for photons emitted from dust and subsequently scattered.

See :ref:`post-processing` for information on how to extract this information
from the output files.

Disabling SEDs or Images
------------------------

When adding a set of binned or peeled images, it is possible to disable the
SED or image part::

    image = m.add_binned_images()  # Images and SEDs
    image = m.add_binned_images(image=False)  # SEDs
    image = m.add_binned_images(sed=False)  # Images

    image = m.add_peeled_images()  # Images and SEDs
    image = m.add_peeled_images(image=False)  # SEDs
    image = m.add_peeled_images(sed=False)  # Images

Example
-------

The following example creates two sets of peeled SEDs/images. The first is used to produce an SED with 250 wavelengths from 0.01 to 5000. microns with uncertainties, and the second is used to produce images at 5 wavelengths between 10 and 100 microns, with image size 100x100 and extending +/-1pc in each direction::

    image1 = m.add_peeled_images(image=False)
    image1.set_wavelength_range(250, 0.01, 5000.)
    image1.set_uncertainties(True)

    image2 = m.add_peeled_images(sed=False)
    image2.set_wavelength_range(5, 10., 100.)
    image2.set_image_size(100, 100)
    image2.set_image_limits(-pc, +pc, -pc, +pc)

