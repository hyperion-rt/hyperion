#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

/* Define constants */
double const LN10 = 2.30258509299404590109;

/* Define docstrings */
static char module_docstring[] = "Fast trapezium integration";
static char integrate_docstring[] = "Integrate in linear space";
static char integrate_loglin_docstring[] = "Integrate in log-linear space";
static char integrate_linlog_docstring[] = "Integrate in linear-log space";
static char integrate_loglog_docstring[] = "Integrate in log-log space";

/* Declare the C functions here. */
static PyObject *_integrate(PyObject *self, PyObject *args);
static PyObject *_integrate_loglin(PyObject *self, PyObject *args);
static PyObject *_integrate_linlog(PyObject *self, PyObject *args);
static PyObject *_integrate_loglog(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
    {"_integrate", _integrate, METH_VARARGS, integrate_docstring},
    {"_integrate_loglin", _integrate_loglin, METH_VARARGS, integrate_loglin_docstring},
    {"_integrate_linlog", _integrate_linlog, METH_VARARGS, integrate_linlog_docstring},
    {"_integrate_loglog", _integrate_loglog, METH_VARARGS, integrate_loglog_docstring},
    {NULL, NULL, 0, NULL}
};

/* This is the function that is called on import. */

#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) void init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(_integrate_core)
{
    PyObject *m;
    MOD_DEF(m, "_integrate_core", module_docstring, module_methods);
    if (m == NULL)
        return MOD_ERROR_VAL;
    import_array();
    return MOD_SUCCESS_VAL(m);
}

/* Do the heavy lifting here */

static PyObject *_integrate(PyObject *self, PyObject *args)
{
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO", &x_obj, &y_obj))
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

    /* Calculate integral */

    double integral = 0.;

    int i, j;

    for (i = 0; i < n - 1; i++)
    {
        j = i + 1;
        if(!npy_isnan(y[i]) && !npy_isnan(y[j])) {
            if(x[j] > x[i]) {
                integral += 0.5 * (x[j] - x[i]) * (y[j] + y[i]);
            } else {
                PyErr_SetString(PyExc_ValueError, "x is not monotonically increasing");
                return NULL;
            }
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", integral);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;
}

static PyObject *_integrate_loglin(PyObject *self, PyObject *args)
{
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO", &x_obj, &y_obj))
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

    /* Calculate integral */

    double integral = 0.;

    int i, j;
    double a, b;

    for (i = 0; i < n - 1; i++)
    {
        j = i + 1;
        if(!npy_isnan(y[i]) && !npy_isnan(y[j])) {
            if(x[j] > x[i]) {
                a = (y[i] - y[j]) / log10(x[i] / x[j]);
                b = y[i] - a * log10(x[i]);
                integral += a * (x[j] * log10(x[j]) - x[i] * log10(x[i])) + (b - a / LN10) * (x[j] - x[i]);
            } else {
                PyErr_SetString(PyExc_ValueError, "x is not monotonically increasing");
                return NULL;
            }
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", integral);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;
}

static PyObject *_integrate_linlog(PyObject *self, PyObject *args)
{
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO", &x_obj, &y_obj))
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

    /* Calculate integral */

    double integral = 0.;

    int i, j;

    for (i = 0; i < n - 1; i++)
    {
        j = i + 1;
        if(!npy_isnan(y[i]) && !npy_isnan(y[j])) {
            if(x[j] > x[i]) {
                if(y[i] == y[j]) {
                    integral += y[i] * (x[j] - x[i]);
                } else {
                    integral += (y[j] - y[i]) * (x[j] - x[i]) / LN10 / log10(y[j] / y[i]);
                }
            } else {
                PyErr_SetString(PyExc_ValueError, "x is not monotonically increasing");
                return NULL;
            }
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", integral);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;
}

static PyObject *_integrate_loglog(PyObject *self, PyObject *args)
{
    PyObject *x_obj, *y_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO", &x_obj, &y_obj))
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

    /* Calculate integral */

    double integral = 0.;

    int i, j;
    double b;

    for (i = 0; i < n - 1; i++)
    {
        j = i + 1;
        if(y[i] > 0. && y[j] > 0. && !npy_isnan(y[i]) && !npy_isnan(y[j])) {
            if(x[j] > x[i]) {
                b = log10(y[i] / y[j]) / log10(x[i] / x[j]);
                if(fabs(b + 1.) < 1.e-10) {
                    integral += x[i] * y[i] * log(x[j] / x[i]);
                } else {
                    integral += y[i] * (x[j] * pow(x[j] / x[i], b) - x[i]) / (b + 1.);
                }
            } else {
                PyErr_SetString(PyExc_ValueError, "x is not monotonically increasing");
                return NULL;
            }
        }
    }

    /* Clean up. */
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", integral);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;
}
