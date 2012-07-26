=================
Writing out files
=================

.. _Numpy: http://numpy.scipy.org/
.. _PyWCS: https://trac6.assembla.com/astrolib
.. _PyFITS: http://www.stsci.edu/resources/software_hardware/pyfits

The output files from the radiative transfer code are in the HDF5 file format,
and can therefore be accessed directly from most programming/scripting
languages. However, in many cases it might be most convenient to write a small
Python script to extract the required information and to write it out to files
that can then be read in to other tools. In this tutorial, we learn how to
write out SEDs and images to ASCII and FITS files respectively.

The tutorial here assumes that you are using the :ref:`tutorial-model`.

Writing files in Python
=======================

Pure Python
-----------

The most basic way to write files in Python is to simply open a file with
write access::

    f = open('file.txt', 'wb')

and to then call the ``write`` method to write to the file::

    f.write("Hello World")

Line returns have to be explicitly included using ``\n``::

    f.write("Line 1\n")
    f.write("line 2\n")

And files should be closed with::

    f.close()

The best way to write out variables with this technique is to use
string formatting which is described in more detail
`here <http://docs.python.org/library/stdtypes.html#string-formatting>`_.
The basic command to format variables into a string is::

    format % variables

where ``format`` is a string containing the format statements and variables is
a tuple of the values, for example::

    >>> print "%s %5.2f %10.4e" % ("name", 3.4, 1.e-10)
    name  3.40 1.0000e-10

We can use this when writing out files, so if we have two lists or arrays of
values ``a`` and ``b`` we can do::

    a = [1,2,3,4,5]
    b = [2,6,4,3,2]

    f = open('file.txt', 'wb')
    for i in range(len(a)):
        f.write("%i %5.2f\n" % (a[i], b[i]))
    f.close()

which will produce a file containing::

    1  2.00
    2  6.00
    3  4.00
    4  3.00
    5  2.00

Numpy
-----

`Numpy`_ provides a function called ``savetxt`` that makes it easy to write out
arrays to files. Given two lists or arrays ``a`` and ``b`` as above, one can
simply do::

   import numpy as np
   a = [1,2,3,4,5]
   b = [2,6,4,3,2]
   np.savetxt('file_numpy.txt', zip(a, b), fmt="%i %5.2f")

which produces exactly the same output as above and avoids the for loop.

Writing out SEDs
================

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

Writing out images
==================

Writing out images to text files does not make much sense, so in this section
we see how to write out images extracted from the radiative transfer code
results to a FITS file, and add WCS information. The first step is to extract
the images from the radiative transfer code. This step is described in detail
in :ref:`post-processing`. Once a 2D image or 3D wavelength cube have been
extracted, we can write them out to a FITS file using
`PyFITS`_::

    import pyfits

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
    pyfits.writeto('image_cube.fits', nufnu.swapaxes(0, 2).swapaxes(1, 2), \
                   clobber=True)

    # We can also just output one of the wavelengths
    pyfits.writeto('image_slice.fits', nufnu[:, :, 0], clobber=True)

Writing out images (with WCS)
=============================

Adding World Coordinate System (WCS) information is easy using
`PyWCS`_::

    import pywcs
    import pyfits

    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc

    m = ModelOutput('tutorial_model.rtout')
    wav, nufnu = m.get_image(group=1, inclination=0, distance=300 * pc)

    # Initialize WCS information
    wcs = pywcs.WCS(naxis=2)

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
    pyfits.writeto('image_slice_wcs.fits', nufnu[:, :, 0], header,
                   clobber=True)
