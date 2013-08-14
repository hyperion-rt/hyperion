import numpy as np

from hyperion.model import ModelOutput
from hyperion.util.constants import pc

m = ModelOutput('class2_sed.rtout')
sed = m.get_sed(inclination=0, aperture=-1, distance=300 * pc)
np.savetxt('sed.txt', list(zip(sed.wav, sed.val)), fmt="%11.4e %11.4e")
