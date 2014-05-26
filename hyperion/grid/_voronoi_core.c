#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#include <math.h>

/* Define docstrings */
static char module_docstring[] = "C implementation of utility functions used in Voronoi grids";
static char simplex3d_volume_docstring[] = "Compute the volume of a 3d simplex";

/* Declare the C functions here. */
static PyObject *_simplex3d_volume(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
    {"_simplex3d_volume", _simplex3d_volume, METH_VARARGS, simplex3d_volume_docstring},
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

MOD_INIT(_voronoi_core)
{
    PyObject *m;
    MOD_DEF(m, "_voronoi_core", module_docstring, module_methods);
    if (m == NULL)
        return MOD_ERROR_VAL;
    import_array();
    return MOD_SUCCESS_VAL(m);
}

static PyObject *_simplex3d_volume(PyObject *self, PyObject *args)
{
    PyObject *arr_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "O", &arr_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *array = PyArray_FROM_OTF(arr_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(array);
        return NULL;
    }

    /* Check the dimensions. */
    int n = (int)PyArray_DIM(array, 0), m = (int)PyArray_DIM(array, 1);

    if (n != 4 || m != 3) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid array dimension(s).");
        Py_DECREF(array);
        return NULL;
    }

    /* Get pointer to the data. */
    double *data = (double*)PyArray_DATA(array);

    /* Build the volume matrix from the input simplex */
    double vmatrix[9];

    int j;
    for (j = 0; j < 3; ++j) {
        vmatrix[j]     = data[3 + j] - data[j];
        vmatrix[j + 3] = data[6 + j] - data[j];
        vmatrix[j + 6] = data[9 + j] - data[j];
    }

    /* Compute the volume. */
    double a,b,c,d,e,f,g,h,i;

    a = vmatrix[0];
    b = vmatrix[1];
    c = vmatrix[2];

    d = vmatrix[3];
    e = vmatrix[4];
    f = vmatrix[5];

    g = vmatrix[6];
    h = vmatrix[7];
    i = vmatrix[8];

    double retval = fabs(a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)) / 6.;

    /* Clean up. */
    Py_DECREF(array);

    /* Build the output value. */
    PyObject *ret = Py_BuildValue("d", retval);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;
}
