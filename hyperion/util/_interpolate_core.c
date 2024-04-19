#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define Py_LIMITED_API 0x030900f0

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

/* Workaround for gcc<10 */
struct _typeobject {int _placeholder;};

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

#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods) \
        static struct PyModuleDef moduledef = { \
        PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
        ob = PyModule_Create(&moduledef);

MOD_INIT(_interpolate_core)
{
    PyObject *m;
    MOD_DEF(m, "_interpolate_core", module_docstring, module_methods);
    if (m == NULL)
        return MOD_ERROR_VAL;
    import_array();
    return MOD_SUCCESS_VAL(m);
}

/* Do the heavy lifting here */

static int binary_locate(double *x, int n, double xval) {

    int imin, imax, imid;

    imin = 0;
    imax = n - 1;

    while(imax > imin + 1) {
        imid = (imin + imax) / 2;
        if(xval > x[imid]) {
            imin = imid;
        } else {
            imax = imid;
        }
    }

    return imax;

}

static PyObject *interp1d_linear_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOd", &x_obj, &y_obj, &xval))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

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

    int ipos = binary_locate(x, n, xval);

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

    PyObject *x_obj, *y_obj, *xval_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOO", &x_obj, &y_obj, &xval_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
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
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2, ipos;

    for(i = 0; i < nval; i++) {
        ipos = binary_locate(x, n, xval[i]);
        if (ipos == 0) {
            yval[i] = y[0];
        } else if (ipos == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos - 1;
            i2 = ipos;
            yval[i] = (xval[i] - x[i1]) / (x[i2] - x[i1]) * (y[i2] - y[i1]) + y[i1];
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);
    Py_DECREF(xval_array);

    return yval_array;

}

static PyObject *interp1d_loglog_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOd", &x_obj, &y_obj, &xval))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

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

    int ipos = binary_locate(x, n, xval);

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

    PyObject *x_obj, *y_obj, *xval_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOO", &x_obj, &y_obj, &xval_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
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
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2, ipos;

    for(i = 0; i < nval; i++) {
        ipos = binary_locate(x, n, xval[i]);
        if (ipos == 0) {
            yval[i] = y[0];
        } else if (ipos == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos - 1;
            i2 = ipos;
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

    return yval_array;

}

static PyObject *interp1d_linlog_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOd", &x_obj, &y_obj, &xval))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

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

    int ipos = binary_locate(x, n, xval);

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

    PyObject *x_obj, *y_obj, *xval_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOO", &x_obj, &y_obj, &xval_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
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
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2, ipos;

    for(i = 0; i < nval; i++) {
        ipos = binary_locate(x, n, xval[i]);
        if (ipos == 0) {
            yval[i] = y[0];
        } else if (ipos == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos - 1;
            i2 = ipos;
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

    return yval_array;

}

static PyObject *interp1d_loglin_scalar(PyObject *self, PyObject *args)
{

    double xval, yval;
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOd", &x_obj, &y_obj, &xval))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

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

    int ipos = binary_locate(x, n, xval);

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

    PyObject *x_obj, *y_obj, *xval_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOO", &x_obj, &y_obj, &xval_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *xval_array = PyArray_FROM_OTF(xval_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (x_array == NULL || y_array == NULL || xval_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(x_array);
        Py_XDECREF(y_array);
        Py_XDECREF(xval_array);
        return NULL;
    }

    /* How many data points are there? */
    int n = (int)PyArray_DIM(x_array, 0);

    /* How many values to interpolate? */
    int nval = (int)PyArray_DIM(xval_array, 0);

    /* Check the dimensions. */
    if (n != (int)PyArray_DIM(y_array, 0)) {
        PyErr_SetString(PyExc_RuntimeError, "Dimension mismatch.");
        Py_DECREF(x_array);
        Py_DECREF(y_array);
        Py_DECREF(xval_array);
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
        Py_XDECREF(yval_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *x = (double*)PyArray_DATA(x_array);
    double *y = (double*)PyArray_DATA(y_array);
    double *xval = (double*)PyArray_DATA(xval_array);
    double *yval = (double*)PyArray_DATA(yval_array);

    /* Interpolate */

    int i, i1, i2, ipos;

    for(i = 0; i < nval; i++) {
        ipos = binary_locate(x, n, xval[i]);
        if (ipos == 0) {
            yval[i] = y[0];
        } else if (ipos == n) {
            yval[i] = y[-1];
        } else {
            i1 = ipos - 1;
            i2 = ipos;
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

    return yval_array;

}
