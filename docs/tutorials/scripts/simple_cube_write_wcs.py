import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

# Retrieve image cube as before
m = ModelOutput('simple_cube.rtout')
image = m.get_image(inclination=0, distance=300 * pc, units='MJy/sr')

# Initialize WCS information
wcs = WCS(naxis=2)

# Use the center of the image as projection center
wcs.wcs.crpix = [image.val.shape[1] / 2. + 0.5,
                 image.val.shape[0] / 2. + 0.5]

# Set the coordinates of the image center
wcs.wcs.crval = [233.4452, 1.2233]

# Set the pixel scale (in deg/pix)
scale = np.degrees(3. * pc / image.val.shape[0] / image.distance)
wcs.wcs.cdelt = [-scale, scale]

# Set the coordinate system
wcs.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']

# And produce a FITS header
header = wcs.to_header()

# We can also just output one of the wavelengths
fits.writeto('simple_cube_slice_wcs.fits', image.val[:, :, 1],
             header=header, clobber=True)
