from astropy.io import fits

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

# Retrieve image cube as before
m = ModelOutput('simple_cube.rtout')
image = m.get_image(inclination=0, distance=300 * pc, units='MJy/sr')

# The image extracted above is a 3D array. We can write it out to FITS.
# We need to swap some of the directions around so as to be able to use
# the ds9 slider to change the wavelength of the image.
fits.writeto('simple_cube.fits', image.val.swapaxes(0, 2).swapaxes(1, 2),
             clobber=True)

# We can also just output one of the wavelengths
fits.writeto('simple_cube_slice.fits', image.val[:, :, 1], clobber=True)
