=================
Writing out files
=================

.. _Numpy: http://numpy.scipy.org/
.. _Astropy: http://www.astropy.org

The output files from the radiative transfer code are in the HDF5 file format,
and can therefore be accessed directly from most programming/scripting
languages. However, in many cases it might be most convenient to write a small
Python script to extract the required information and to write it out to files
that can then be read in to other tools. In this tutorial, we learn how to
write out SEDs and images to ASCII and FITS files respectively.

The tutorial here assumes that you are using the :ref:`tutorial-model`.

.. note:: If you are not familiar with writing files from Python, you can first 
          take a look at the :doc:`python_writing` tutorial.


SEDs
====

We are now ready to write out SEDs to ASCII files. The first step is to
extract the SED from the output file from the radiation transfer code. This
step is described in detail in :ref:`post-processing`. Combining this with what
we learned above about writing files, we can write scripts that will fetch
SEDs and write them out to disk. For example, if we want to write out the SED
for the first inclination and the largest aperture, we can do::

    import numpy as np

    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc

    # Open the model - we specify the name without the .rtout extension
    m = ModelOutput('tutorial_model.rtout')

    # Extract the SED for the smallest inclination and largest aperture, and
    # scale to 300pc. In Python, negative indices can be used for lists and
    # arrays, and indicate the position from the end. So to get the SED in the
    # largest aperture, we set aperture=-1.
    wav, nufnu = m.get_sed(inclination=0, aperture=-1, distance=300 * pc)

    # Write out the SED to file
    np.savetxt('sed.txt', zip(wav, nufnu), fmt="%11.4e %11.4e")

This script produces a file that looks like::

    4.8705e+03  7.5202e-13
    4.6214e+03  9.6952e-13
    4.3851e+03  1.2493e-12
    4.1609e+03  1.6091e-12
    3.9481e+03  2.0714e-12
    3.7462e+03  2.6649e-12
    3.5546e+03  3.4265e-12
    3.3729e+03  4.4029e-12
    ...

Images
======

Writing out images to text files does not make much sense, so in this section
we see how to write out images extracted from the radiative transfer code
results to a FITS file, and add WCS information. The first step is to extract
the images from the radiative transfer code. This step is described in detail
in :ref:`post-processing`. Once a 2D image or 3D wavelength cube have been
extracted, we can write them out to a FITS file using
`Astropy`_::

    from astropy.io import fits

    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc

    # Open the model - we specify the name without the .rtout extension
    m = ModelOutput('tutorial_model.rtout')

    # Extract the image for the first inclination, and scale to 300pc. We
    # have to specify group=1 as there is no image in group 0
    wav, nufnu = m.get_image(group=1, inclination=0, distance=300 * pc)

    # The image extracted above is a 3D array. We can write it out to FITS.
    # We need to swap some of the directions around so as to be able to use
    # the ds9 slider to change the wavelength of the image.
    fits.writeto('image_cube.fits', nufnu.swapaxes(0, 2).swapaxes(1, 2), \
                   clobber=True)

    # We can also just output one of the wavelengths
    fits.writeto('image_slice.fits', nufnu[:, :, 0], clobber=True)

Images (with WCS)
=================

Adding World Coordinate System (WCS) information is easy using
`Astropy`_::

    from astropy.io import fits
    from astropy.wcs import WCS

    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc

    m = ModelOutput('tutorial_model.rtout')
    wav, nufnu = m.get_image(group=1, inclination=0, distance=300 * pc)

    # Initialize WCS information
    wcs = WCS(naxis=2)

    # Use the center of the image as projection center
    wcs.wcs.crpix = [nufnu.shape[2] / 2. + 0.5,
                     nufnu.shape[1] / 2. + 0.5]

    # Set the coordinates of the image center
    wcs.wcs.crval = [233.4452, 1.2233]

    # Set the pixel scale (in deg/pix)
    wcs.wcs.cdelt = [1./3600., 1./3600.]

    # Set the coordinate system
    wcs.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']

    # And produce a FITS header
    header = wcs.to_header()

    # Write out to a file including the new header
    fits.writeto('image_slice_wcs.fits', nufnu[:, :, 0], header,
                 clobber=True)
