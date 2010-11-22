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
        return self.value * rstar \
               * (1. - (1. - 2. * (self.temperature / tstar) ** 4. \
               * dust.kappa_planck_temperature(self.temperature) \
               / dust.kappa_planck_spectrum(nu, fnu)) ** 2.) ** -0.5
