import numpy as np

from ..functions import B_nu, dB_nu_dT
from ..integrate import integrate_loglog
from ..constants import sigma, k, c


def test_b_nu():

    nu = np.logspace(-20, 20., 10000)

    for T in [10,100,1000,10000]:

        # Compute planck function
        b = B_nu(nu, T)

        # Check that the intergral is correct
        total = integrate_loglog(nu, b)
        np.testing.assert_allclose(total, sigma * T ** 4 / np.pi, rtol=1e-4)

        # Check that we reach the rayleigh-jeans limit at low frequencies
        rj = 2. * nu **2 * k * T / c**2
        np.testing.assert_allclose(b[nu < 1e-10], rj[nu < 1e-10], rtol=1.e-8)


def test_db_nu_dt():

    nu = np.logspace(-20, 20., 10000)

    for T in [10,100,1000,10000]:

        # Compute exact planck function derivative
        db = dB_nu_dT(nu, T)

        # Compute numerical planck function derivative
        dT = T / 1e6
        b1 = B_nu(nu, T - dT)
        b2 = B_nu(nu, T + dT)
        db_num = 0.5 * (b2 - b1) / dT

        # Check that the two are the same
        np.testing.assert_allclose(db, db_num, rtol=1.e-2)
