import numpy as np
from PIL import Image

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

m = ModelOutput('simple_cube.rtout')
image = m.get_image(inclination=0, distance=300 * pc, units='MJy/sr')

# Extract the slices we want to use for red, green, and blue
r = image.val[:, :, 17]
g = image.val[:, :, 18]
b = image.val[:, :, 19]

# Now we need to rescale the values we want to the range 0 to 255, clip values
# outside the range, and convert to unsigned 8-bit integers. We also use a sqrt
# stretch (hence the ** 0.5)

r = np.clip((r / 0.5) ** 0.5 * 255., 0., 255.)
r = np.array(r, dtype=np.uint8)

g = np.clip((g / 2) ** 0.5 * 255., 0., 255.)
g = np.array(g, dtype=np.uint8)

b = np.clip((b / 4.) ** 0.5 * 255., 0., 255.)
b = np.array(b, dtype=np.uint8)

# We now convert to image objects
image_r = Image.fromarray(r)
image_g = Image.fromarray(g)
image_b = Image.fromarray(b)

# And finally merge into a single 3-color image
img = Image.merge("RGB", (image_r, image_g, image_b))

# By default, the image will be flipped, so we need to fix this
img = img.transpose(Image.FLIP_TOP_BOTTOM)

img.save('simple_cube_rgb.png')
