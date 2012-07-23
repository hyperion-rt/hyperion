#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

/* Define docstrings */
static char module_docstring[] = "Fast trapezium integration";
static char interp1d_linear_scalar_docstring[] = "Interpolate in linear-linear space";
static char interp1d_linear_array_docstring[] = "Interpolate in linear-linear space (array)";
static char interp1d_loglog_scalar_docstring[] = "Interpolate in log-log space";
static char interp1d_loglog_array_docstring[] = "Interpolate in log-log space (array)";
static char interp1d_linlog_scalar_docstring[] = "Interpolate in linear-log space";
static char interp1d_linlog_array_docstring[] = "Interpolate in linear-log space (array)";
static char interp1d_loglin_scalar_docstring[] = "Interpolate in log-linear space";
static char interp1d_loglin_array_docstring[] = "Interpolate in log-linear space (array)";

/* Declare the C functions here. */
static PyObject *interp1d_linear_scalar(PyObject *self, PyObject *args);
static PyObject *interp1d_linear_array(PyObject *self, PyObject *args);
static PyObject *interp1d_loglog_scalar(PyObject *self, PyObject *args);
static PyObject *interp1d_loglog_array(PyObject *self, PyObject *args);
static PyObject *interp1d_linlog_scalar(PyObject *self, PyObject *args);
static PyObject *interp1d_linlog_array(PyObject *self, PyObject *args);
static PyObject *interp1d_loglin_scalar(PyObject *self, PyObject *args);
static PyObject *interp1d_loglin_array(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
    {"interp1d_linear_scalar", interp1d_linear_scalar, METH_VARARGS, interp1d_linear_scalar_docstring},
    {"interp1d_linear_array", interp1d_linear_array, METH_VARARGS, interp1d_linear_array_docstring},
    {"interp1d_loglog_scalar", interp1d_loglog_scalar, METH_VARARGS, interp1d_loglog_scalar_docstring},
    {"interp1d_loglog_array", interp1d_loglog_array, METH_VARARGS, interp1d_loglog_array_docstring},
    {"interp1d_linlog_scalar", interp1d_linlog_scalar, METH_VARARGS, interp1d_linlog_scalar_docstring},
    {"interp1d_linlog_array", interp1d_linlog_array, METH_VARARGS, interp1d_linlog_array_docstring},
    {"interp1d_loglin_scalar", interp1d_loglin_scalar, METH_VARARGS, interp1d_loglin_scalar_docstring},
    {"interp1d_loglin_array", interp1d_loglin_array, METH_VARARGS, interp1d_loglin_array_docstring},
    {NULL, NULL, 0, NULL}
};

/* This is the function that is called on import. */

PyMODINIT_FUNC init_interpolate_core(void)
{
    /* Initialize the module with a docstring. */
    PyObject *m = Py_InitModule3("_interpolate_core", module_methods, module_docstring);
    if (m == NULL)
        return;

    /* Load all of the `numpy` functionality. */
    import_array();
}

/* Do the heavy lifting here */

static PyObject *interp1d_linear_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    int ipos;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOdi", &x_obj, &y_obj, &xval, &ipos))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);

    /* Interpolate */

    int i1, i2;

    if (ipos == 0) {
        yval = y[0];
    } else if (ipos == n) {
        yval = y[-1];
    } else {
        i1 = ipos - 1;
        i2 = ipos;
        yval = (xval - x[i1]) / (x[i2] - x[i1]) * (y[i2] - y[i1]) + y[i1];
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", yval);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;

}

static PyObject *interp1d_linear_array(PyObject *self, PyObject *args)
{

    PyObject *x_obj, *y_obj, *xval_obj, *ipos_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOO", &x_obj, &y_obj, &xval_obj, &ipos_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ipos_array = PyArray_FROM_OTF(ipos_obj, NPY_INT64, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL || ipos_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        Py_XDECREF(ipos_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0) || nval != (int)PyArray_DIM(ipos_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_XDECREF(xval_array);
        Py_XDECREF(ipos_array);
        return NULL;
    }

    /* Build the output array */
    npy_intp dims[1];
    dims[0] = nval;
    PyObject *yval_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (yval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
        Py_DECREF(ipos_array);
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    long *ipos = (long*)PyArray_DATA(ipos_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2;

    for(i = 0; i < nval; i++) {
        if (ipos[i] == 0) {
            yval[i] = y[0];
        } else if (ipos[i] == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos[i] - 1;
            i2 = ipos[i];
            yval[i] = (xval[i] - x[i1]) / (x[i2] - x[i1]) * (y[i2] - y[i1]) + y[i1];
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(xval_array);
    Py_DECREF(ipos_array);

    return yval_array;

}

static PyObject *interp1d_loglog_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    int ipos;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOdi", &x_obj, &y_obj, &xval, &ipos))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);

    /* Interpolate */

    int i1, i2;

    if (ipos == 0) {
        yval = y[0];
    } else if (ipos == n) {
        yval = y[-1];
    } else {
        i1 = ipos - 1;
        i2 = ipos;
        if (y[i1] > 0. && y[i2] > 0.) {
            yval = pow(10., ((log10(xval) - log10(x[i1])) \
                           / (log10(x[i2]) - log10(x[i1])) \
                           * (log10(y[i2]) - log10(y[i1])) \
                           + log10(y[i1])));
        } else {
            yval = 0.;
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", yval);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;

}

static PyObject *interp1d_loglog_array(PyObject *self, PyObject *args)
{

    PyObject *x_obj, *y_obj, *xval_obj, *ipos_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOO", &x_obj, &y_obj, &xval_obj, &ipos_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ipos_array = PyArray_FROM_OTF(ipos_obj, NPY_INT64, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL || ipos_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        Py_XDECREF(ipos_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0) || nval != (int)PyArray_DIM(ipos_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
        Py_DECREF(ipos_array);
        return NULL;
    }

    /* Build the output array */
    npy_intp dims[1];
    dims[0] = nval;
    PyObject *yval_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (yval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
        Py_DECREF(ipos_array);
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    long *ipos = (long*)PyArray_DATA(ipos_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2;

    for(i = 0; i < nval; i++) {
        if (ipos[i] == 0) {
            yval[i] = y[0];
        } else if (ipos[i] == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos[i] - 1;
            i2 = ipos[i];
            if (y[i1] > 0. && y[i2] > 0.) {
                yval[i] = pow(10., ((log10(xval[i]) - log10(x[i1])) \
                                  / (log10(x[i2]) - log10(x[i1])) \
                                  * (log10(y[i2]) - log10(y[i1])) \
                                  + log10(y[i1])));
            } else {
                yval[i] = 0.;
            }
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(xval_array);
    Py_DECREF(ipos_array);

    return yval_array;

}

static PyObject *interp1d_linlog_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    int ipos;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOdi", &x_obj, &y_obj, &xval, &ipos))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);

    /* Interpolate */

    int i1, i2;

    if (ipos == 0) {
        yval = y[0];
    } else if (ipos == n) {
        yval = y[-1];
    } else {
        i1 = ipos - 1;
        i2 = ipos;
        if (y[i1] > 0. && y[i2] > 0.) {
            yval = pow(10., (xval - x[i1]) \
                           / (x[i2] - x[i1]) \
                           * (log10(y[i2]) - log10(y[i1])) \
                           + log10(y[i1]));
        } else {
            yval = 0.;
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", yval);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;

}

static PyObject *interp1d_linlog_array(PyObject *self, PyObject *args)
{

    PyObject *x_obj, *y_obj, *xval_obj, *ipos_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOO", &x_obj, &y_obj, &xval_obj, &ipos_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ipos_array = PyArray_FROM_OTF(ipos_obj, NPY_INT64, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL || ipos_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        Py_XDECREF(ipos_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0) || nval != (int)PyArray_DIM(ipos_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
        Py_DECREF(ipos_array);
        return NULL;
    }

    /* Build the output array */
    npy_intp dims[1];
    dims[0] = nval;
    PyObject *yval_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (yval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
        Py_DECREF(ipos_array);
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    long *ipos = (long*)PyArray_DATA(ipos_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2;

    for(i = 0; i < nval; i++) {
        if (ipos[i] == 0) {
            yval[i] = y[0];
        } else if (ipos[i] == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos[i] - 1;
            i2 = ipos[i];
            if (y[i1] > 0. && y[i2] > 0.) {
                yval[i] = pow(10., (xval[i] - x[i1]) \
                                  / (x[i2] - x[i1]) \
                                  * (log10(y[i2]) - log10(y[i1])) \
                                  + log10(y[i1]));
            } else {
                yval[i] = 0.;
            }
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(xval_array);
    Py_DECREF(ipos_array);

    return yval_array;

}

static PyObject *interp1d_loglin_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    int ipos;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOdi", &x_obj, &y_obj, &xval, &ipos))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);

    /* Interpolate */

    int i1, i2;

    if (ipos == 0) {
        yval = y[0];
    } else if (ipos == n) {
        yval = y[-1];
    } else {
        i1 = ipos - 1;
        i2 = ipos;
        yval =  ((log10(xval) - log10(x[i1])) \
               / (log10(x[i2]) - log10(x[i1])) \
               * (y[i2] - y[i1]) \
               + y[i1]);
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", yval);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;

}

static PyObject *interp1d_loglin_array(PyObject *self, PyObject *args)
{

    PyObject *x_obj, *y_obj, *xval_obj, *ipos_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOO", &x_obj, &y_obj, &xval_obj, &ipos_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ipos_array = PyArray_FROM_OTF(ipos_obj, NPY_INT64, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL || ipos_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        Py_XDECREF(ipos_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0) || nval != (int)PyArray_DIM(ipos_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
        Py_DECREF(ipos_array);
        return NULL;
    }

    /* Build the output array */
    npy_intp dims[1];
    dims[0] = nval;
    PyObject *yval_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (yval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
        Py_DECREF(ipos_array);
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    long *ipos = (long*)PyArray_DATA(ipos_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2;

    for(i = 0; i < nval; i++) {
        if (ipos[i] == 0) {
            yval[i] = y[0];
        } else if (ipos[i] == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos[i] - 1;
            i2 = ipos[i];
            yval[i] =  ((log10(xval[i]) - log10(x[i1])) \
                      / (log10(x[i2]) - log10(x[i1])) \
                      * (y[i2] - y[i1]) \
                      + y[i1]);
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(xval_array);
    Py_DECREF(ipos_array);

    return yval_array;

}