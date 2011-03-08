import numpy as np


class OptThinRadius(object):

    def __init__(self, temperature, value=1.):
        self.temperature = temperature
        self.value = value

    def __mul__(self, value):
        return OptThinRadius(self.temperature, value=self.value * value)

    def __rmul__(self, value):
        return OptThinRadius(self.temperature, value=self.value * value)

    def __str__(self):
        return "%g times the dust sublimation radius" % self.n

    def evaluate(self, star, dust):
        rstar = star.radius
        tstar = star.effective_temperature()
        nu, fnu = star.total_spectrum()
        x = (self.temperature / tstar) ** 4. \
            * dust.kappa_planck_temperature(self.temperature) \
            / dust.kappa_planck_spectrum(nu, fnu)
        if x < 0.001:
            r = self.value * rstar / 2. / np.sqrt(x)
        else:
            r = self.value * rstar / np.sqrt(1. - (1. - 2. * x) ** 2.)
        if np.isnan(r):
            raise Exception("Optically thin radius is NaN")
        return r
