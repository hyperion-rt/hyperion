from __future__ import print_function, division

import numpy as np


class OptThinRadius(object):

    def __init__(self, temperature, value=1., min=0.):
        self.temperature = temperature
        self.value = value
        self.min = 0.

    def __mul__(self, value):
        return OptThinRadius(self.temperature, value=self.value * value)

    def __rmul__(self, value):
        return OptThinRadius(self.temperature, value=self.value * value)

    def __str__(self):
        return "%g times the radius at which the optically thin temperature would be %gK" % (self.value, self.temperature)

    def evaluate(self, star, dust):
        rstar = star.radius
        tstar = star.effective_temperature()
        nu_min, nu_max = dust.optical_properties.nu[0], \
                         dust.optical_properties.nu[-1]
        nu, fnu = star.total_spectrum(bnu_range=(nu_min, nu_max))
        x = (self.temperature / tstar) ** 4. \
            * dust.optical_properties.kappa_planck_temperature(self.temperature) \
            / dust.optical_properties.kappa_planck_spectrum(nu, fnu)
        if x < 0.001:
            r = self.value * rstar / 2. / np.sqrt(x)
        else:
            r = self.value * rstar / np.sqrt(1. - (1. - 2. * x) ** 2.)
        if np.isnan(r):
            raise Exception("Optically thin radius is NaN")
        return max(r, self.min)
