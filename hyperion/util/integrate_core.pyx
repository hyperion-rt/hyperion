from __future__ import division
import numpy as np
cimport numpy as np

from libc.math cimport log, log10, pow

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "numpy/npy_math.h":
    bint npy_isnan(double x)

cimport cython

LN10 = log(10.)


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def integrate(np.ndarray[DTYPE_t, ndim=1] x,
              np.ndarray[DTYPE_t, ndim=1] y):

    assert x.dtype == DTYPE and y.dtype == DTYPE

    cdef int n = x.shape[0]
    cdef DTYPE_t integral
    cdef unsigned int i, j

    integral = 0.
    for i in range(n - 1):
        j = i + 1
        if not npy_isnan(y[i]) and not npy_isnan(y[j]):
            if x[j] > x[i]:
                integral += 0.5 * (x[j] - x[i]) * (y[j] + y[i])
            elif x[j] < x[i]:
                raise ValueError('x is not monotonically increasing')

    return integral


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def integrate_loglin(np.ndarray[DTYPE_t, ndim=1] x,
                     np.ndarray[DTYPE_t, ndim=1] y):

    assert x.dtype == DTYPE and y.dtype == DTYPE

    cdef int n = x.shape[0]
    cdef DTYPE_t integral, a, b
    cdef unsigned int i, j

    integral = 0.
    for i in range(n - 1):
        j = i + 1
        if not npy_isnan(y[i]) and not npy_isnan(y[j]):
            if x[j] > x[i]:
                a = (y[i] - y[j]) / log10(x[i] / x[j])
                b = y[i] - a * log10(x[i])
                integral += a * (x[j] * log10(x[j]) - x[i] * log10(x[i])) + (b - a / LN10) * (x[j] - x[i])
            elif x[j] < x[i]:
                raise ValueError('x is not monotonically increasing')

    return integral


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def integrate_linlog(np.ndarray[DTYPE_t, ndim=1] x,
                     np.ndarray[DTYPE_t, ndim=1] y):

    assert x.dtype == DTYPE and y.dtype == DTYPE

    cdef int n = x.shape[0]
    cdef DTYPE_t integral
    cdef unsigned int i, j

    integral = 0.
    for i in range(n - 1):
        j = i + 1
        if y[i] > 0. and y[j] > 0. and not npy_isnan(y[i]) and not npy_isnan(y[j]):
            if x[j] > x[i]:
                if y[i] == y[j]:
                    integral += y[i] * (x[j] - x[i])
                else:
                    integral += (y[j] - y[i]) * (x[j] - x[i]) / LN10 / log10(y[j] / y[i])
            elif x[j] < x[i]:
                raise ValueError('x is not monotonically increasing')

    return integral


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def integrate_loglog(np.ndarray[DTYPE_t, ndim=1] x,
                     np.ndarray[DTYPE_t, ndim=1] y):

    assert x.dtype == DTYPE and y.dtype == DTYPE

    cdef int n = x.shape[0]
    cdef DTYPE_t integral, b
    cdef unsigned int i, j

    integral = 0.
    for i in range(n - 1):
        j = i + 1
        if y[i] > 0. and y[j] > 0. and not npy_isnan(y[i]) and not npy_isnan(y[j]):
            if x[j] > x[i]:
                b = log10(y[i] / y[j]) / log10(x[i] / x[j])
                if abs(b + 1.) < 1.e-10:
                    integral += x[i] * y[i] * log(x[j] / x[i])
                else:
                    integral += y[i] * (x[j] * pow(x[j] / x[i], b) - x[i]) / (b + 1.)
            elif x[j] < x[i]:
                raise ValueError('x is not monotonically increasing')

    return integral
