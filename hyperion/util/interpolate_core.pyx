from __future__ import division
import numpy as np
cimport numpy as np

DTYPE_F = np.float
DTYPE_I = np.int
ctypedef np.float_t DTYPE_F_t
ctypedef np.int_t DTYPE_I_t

cdef extern from "numpy/npy_math.h":
    bint npy_isnan(double x)
    double npy_log(double x)
    double npy_log10(double x)
    double npy_fabs(double x)
    double npy_pow(double x, double y)

cimport cython


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_linear_scalar(np.ndarray[DTYPE_F_t, ndim=1] x,
                           np.ndarray[DTYPE_F_t, ndim=1] y,
                           double xval,
                           int ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F

    cdef int n = x.shape[0]
    cdef double yval
    cdef unsigned int i1, i2

    if ipos == 0:
        yval = y[0]
    elif ipos == n:
        yval = y[-1]
    else:
        i1 = ipos - 1
        i2 = ipos
        yval = (xval - x[i1]) / (x[i2] - x[i1]) * (y[i2] - y[i1]) + y[i1]

    return yval


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_linear_array(np.ndarray[DTYPE_F_t, ndim=1] x,
                          np.ndarray[DTYPE_F_t, ndim=1] y,
                          np.ndarray[DTYPE_F_t, ndim=1] xval,
                          np.ndarray[DTYPE_I_t, ndim=1] ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F and xval.dtype == DTYPE_F and ipos.dtype == DTYPE_I

    cdef int n = x.shape[0]
    cdef int nval = xval.shape[0]
    cdef np.ndarray[DTYPE_F_t, ndim=1] yval = np.zeros([nval], dtype=DTYPE_F)
    cdef unsigned int i, i1, i2

    for i in range(nval):
        if ipos[i] == 0:
            yval[i] = y[0]
        elif ipos[i] == n:
            yval[i] = y[-1]
        else:
            i1 = ipos[i] - 1
            i2 = ipos[i]
            yval[i] = (xval[i] - x[i1]) / (x[i2] - x[i1]) * (y[i2] - y[i1]) + y[i1]

    return yval


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_loglog_scalar(np.ndarray[DTYPE_F_t, ndim=1] x,
                           np.ndarray[DTYPE_F_t, ndim=1] y,
                           double xval,
                           int ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F

    cdef int n = x.shape[0]
    cdef double yval
    cdef unsigned int i1, i2

    if ipos == 0:
        yval = y[0]
    elif ipos == n:
        yval = y[-1]
    else:
        i1 = ipos - 1
        i2 = ipos
        if y[i1] > 0. and y[i2] > 0.:
            yval = 10. ** ((npy_log10(xval) - npy_log10(x[i1])) \
                         / (npy_log10(x[i2]) - npy_log10(x[i1])) \
                         * (npy_log10(y[i2]) - npy_log10(y[i1])) \
                         + npy_log10(y[i1]))
        else:
            yval = 0.

    return yval


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_loglog_array(np.ndarray[DTYPE_F_t, ndim=1] x,
                          np.ndarray[DTYPE_F_t, ndim=1] y,
                          np.ndarray[DTYPE_F_t, ndim=1] xval,
                          np.ndarray[DTYPE_I_t, ndim=1] ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F and xval.dtype == DTYPE_F and ipos.dtype == DTYPE_I

    cdef int n = x.shape[0]
    cdef int nval = xval.shape[0]
    cdef np.ndarray[DTYPE_F_t, ndim=1] yval = np.zeros([nval], dtype=DTYPE_F)
    cdef unsigned int i, i1, i2

    for i in range(nval):
        if ipos[i] == 0:
            yval[i] = y[0]
        elif ipos[i] == n:
            yval[i] = y[-1]
        else:
            i1 = ipos[i] - 1
            i2 = ipos[i]
            if y[i1] > 0. and y[i2] > 0.:
                yval[i] = 10. ** ((npy_log10(xval[i]) - npy_log10(x[i1])) \
                                / (npy_log10(x[i2]) - npy_log10(x[i1])) \
                                * (npy_log10(y[i2]) - npy_log10(y[i1])) \
                                + npy_log10(y[i1]))
            else:
                yval[i] = 0.

    return yval
    
@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_linlog_scalar(np.ndarray[DTYPE_F_t, ndim=1] x,
                           np.ndarray[DTYPE_F_t, ndim=1] y,
                           double xval,
                           int ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F

    cdef int n = x.shape[0]
    cdef double yval
    cdef unsigned int i1, i2

    if ipos == 0:
        yval = y[0]
    elif ipos == n:
        yval = y[-1]
    else:
        i1 = ipos - 1
        i2 = ipos
        if y[i1] > 0. and y[i2] > 0.:
            yval = 10. ** ((xval - x[i1]) \
                         / (x[i2] - x[i1]) \
                         * (npy_log10(y[i2]) - npy_log10(y[i1])) \
                         + npy_log10(y[i1]))
        else:
            yval = 0.

    return yval


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_linlog_array(np.ndarray[DTYPE_F_t, ndim=1] x,
                          np.ndarray[DTYPE_F_t, ndim=1] y,
                          np.ndarray[DTYPE_F_t, ndim=1] xval,
                          np.ndarray[DTYPE_I_t, ndim=1] ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F and xval.dtype == DTYPE_F and ipos.dtype == DTYPE_I

    cdef int n = x.shape[0]
    cdef int nval = xval.shape[0]
    cdef np.ndarray[DTYPE_F_t, ndim=1] yval = np.zeros([nval], dtype=DTYPE_F)
    cdef unsigned int i, i1, i2

    for i in range(nval):
        if ipos[i] == 0:
            yval[i] = y[0]
        elif ipos[i] == n:
            yval[i] = y[-1]
        else:
            i1 = ipos[i] - 1
            i2 = ipos[i]
            if y[i1] > 0. and y[i2] > 0.:
                yval[i] = 10. ** ((xval[i] - x[i1]) \
                                / (x[i2] - x[i1]) \
                                * (npy_log10(y[i2]) - npy_log10(y[i1])) \
                                + npy_log10(y[i1]))
            else:
                yval[i] = 0.

    return yval
    
@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_loglin_scalar(np.ndarray[DTYPE_F_t, ndim=1] x,
                           np.ndarray[DTYPE_F_t, ndim=1] y,
                           double xval,
                           int ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F

    cdef int n = x.shape[0]
    cdef double yval
    cdef unsigned int i1, i2

    if ipos == 0:
        yval = y[0]
    elif ipos == n:
        yval = y[-1]
    else:
        i1 = ipos - 1
        i2 = ipos
        yval =  ((npy_log10(xval) - npy_log10(x[i1])) \
               / (npy_log10(x[i2]) - npy_log10(x[i1])) \
               * (y[i2] - y[i1]) \
               + y[i1])

    return yval


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def interp1d_loglin_array(np.ndarray[DTYPE_F_t, ndim=1] x,
                          np.ndarray[DTYPE_F_t, ndim=1] y,
                          np.ndarray[DTYPE_F_t, ndim=1] xval,
                          np.ndarray[DTYPE_I_t, ndim=1] ipos):

    assert x.dtype == DTYPE_F and y.dtype == DTYPE_F and xval.dtype == DTYPE_F and ipos.dtype == DTYPE_I

    cdef int n = x.shape[0]
    cdef int nval = xval.shape[0]
    cdef np.ndarray[DTYPE_F_t, ndim=1] yval = np.zeros([nval], dtype=DTYPE_F)
    cdef unsigned int i, i1, i2

    for i in range(nval):
        if ipos[i] == 0:
            yval[i] = y[0]
        elif ipos[i] == n:
            yval[i] = y[-1]
        else:
            i1 = ipos[i] - 1
            i2 = ipos[i]
            yval[i] =  ((npy_log10(xval[i]) - npy_log10(x[i1])) \
                      / (npy_log10(x[i2]) - npy_log10(x[i1])) \
                      * (y[i2] - y[i1]) \
                      + y[i1])

    return yval
