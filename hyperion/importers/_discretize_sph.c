#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include <math.h>

/* Define docstrings */
static char module_docstring[] = "Helpers for discretizing SPH particles";
static char discretize_sph_docstring[] = "Discretize SPH particles onto given cells";

/* Declare the C functions here. */
static PyObject *_discretize_sph_func(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
    {"_discretize_sph_func", _discretize_sph_func, METH_VARARGS, discretize_sph_docstring},
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

MOD_INIT(_discretize_sph)
{
    PyObject *m;
    MOD_DEF(m, "_discretize_sph", module_docstring, module_methods);
    if (m == NULL)
        return MOD_ERROR_VAL;
    import_array();
    return MOD_SUCCESS_VAL(m);
}

/* Do the heavy lifting here */

static PyObject *_discretize_sph_func(PyObject *self, PyObject *args)
{
        
    PyObject *xmin_obj, *xmax_obj;
    PyObject *ymin_obj, *ymax_obj;
    PyObject *zmin_obj, *zmax_obj;

    PyObject *mu_x_obj, *sigma_x_obj;
    PyObject *mu_y_obj, *sigma_y_obj;
    PyObject *mu_z_obj, *sigma_z_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOO", &xmin_obj, &xmax_obj,
                                                &ymin_obj, &ymax_obj,
                                                &zmin_obj, &zmax_obj,
                                                &mu_x_obj, &mu_y_obj, &mu_z_obj,
                                                &sigma_x_obj, &sigma_y_obj, &sigma_z_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *xmin_array = PyArray_FROM_OTF(xmin_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *xmax_array = PyArray_FROM_OTF(xmax_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ymin_array = PyArray_FROM_OTF(ymin_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ymax_array = PyArray_FROM_OTF(ymax_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *zmin_array = PyArray_FROM_OTF(zmin_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *zmax_array = PyArray_FROM_OTF(zmax_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *mu_x_array = PyArray_FROM_OTF(mu_x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *mu_y_array = PyArray_FROM_OTF(mu_y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *mu_z_array = PyArray_FROM_OTF(mu_z_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *sigma_x_array = PyArray_FROM_OTF(sigma_x_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *sigma_y_array = PyArray_FROM_OTF(sigma_y_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *sigma_z_array = PyArray_FROM_OTF(sigma_z_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (xmin_array == NULL || xmax_array == NULL ||
        ymin_array == NULL || ymax_array == NULL ||
        zmin_array == NULL || zmax_array == NULL ||
        mu_x_array == NULL || mu_y_array == NULL || mu_z_array == NULL ||
        sigma_x_array == NULL || sigma_y_array == NULL || sigma_z_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(xmin_array);
        Py_XDECREF(xmax_array);
        Py_XDECREF(ymin_array);
        Py_XDECREF(ymax_array);
        Py_XDECREF(zmin_array);
        Py_XDECREF(zmax_array);
        Py_XDECREF(mu_x_array);
        Py_XDECREF(mu_y_array);
        Py_XDECREF(mu_z_array);
        Py_XDECREF(sigma_x_array);
        Py_XDECREF(sigma_y_array);
        Py_XDECREF(sigma_z_array);
        return NULL;
    }

    /* How many cells are there? */
    int ncells = (int)PyArray_DIM(xmin_array, 0);

    /* How many particles are there? */
    int nsph = (int)PyArray_DIM(mu_x_array, 0);

    /* Check the dimensions. */
    if (ncells != (int)PyArray_DIM(xmax_array, 0) ||
        ncells != (int)PyArray_DIM(ymin_array, 0) ||
        ncells != (int)PyArray_DIM(ymax_array, 0) ||
        ncells != (int)PyArray_DIM(zmin_array, 0) ||
        ncells != (int)PyArray_DIM(zmax_array, 0) ||
        nsph != (int)PyArray_DIM(mu_y_array, 0) ||
        nsph != (int)PyArray_DIM(mu_z_array, 0) ||
        nsph != (int)PyArray_DIM(sigma_x_array, 0) ||
        nsph != (int)PyArray_DIM(sigma_y_array, 0) ||
        nsph != (int)PyArray_DIM(sigma_z_array, 0)
    ) {
        PyErr_SetString(PyExc_RuntimeError, "xmax array dimension mismatch.");
        Py_XDECREF(xmin_array);
        Py_XDECREF(xmax_array);
        Py_XDECREF(ymin_array);
        Py_XDECREF(ymax_array);
        Py_XDECREF(zmin_array);
        Py_XDECREF(zmax_array);
        Py_XDECREF(mu_x_array);
        Py_XDECREF(mu_y_array);
        Py_XDECREF(mu_z_array);
        Py_XDECREF(sigma_x_array);
        Py_XDECREF(sigma_y_array);
        Py_XDECREF(sigma_z_array);
        return NULL;
    }
    
    /* Build the output array */
    npy_intp dims[1];
    dims[0] = ncells;
    PyObject *total_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (total_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_XDECREF(xmin_array);
        Py_XDECREF(xmax_array);
        Py_XDECREF(ymin_array);
        Py_XDECREF(ymax_array);
        Py_XDECREF(zmin_array);
        Py_XDECREF(zmax_array);
        Py_XDECREF(mu_x_array);
        Py_XDECREF(mu_y_array);
        Py_XDECREF(mu_z_array);
        Py_XDECREF(sigma_x_array);
        Py_XDECREF(sigma_y_array);
        Py_XDECREF(sigma_z_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *xmin = (double*)PyArray_DATA(xmin_array);
    double *xmax = (double*)PyArray_DATA(xmax_array);
    double *ymin = (double*)PyArray_DATA(ymin_array);
    double *ymax = (double*)PyArray_DATA(ymax_array);
    double *zmin = (double*)PyArray_DATA(zmin_array);
    double *zmax = (double*)PyArray_DATA(zmax_array);
    double *mu_x = (double*)PyArray_DATA(mu_x_array);
    double *mu_y = (double*)PyArray_DATA(mu_y_array);
    double *mu_z = (double*)PyArray_DATA(mu_z_array);
    double *sigma_x = (double*)PyArray_DATA(sigma_x_array);
    double *sigma_y = (double*)PyArray_DATA(sigma_y_array);
    double *sigma_z = (double*)PyArray_DATA(sigma_z_array);
    double *total = (double*)PyArray_DATA(total_array);
    
    /* Calculate total */

    int i, j;

    /* Loop over all cells */
    for (i = 0; i < ncells - 1; i++)
    {
        
        /* Loop over all SPH particles */

        for (j = 0; j < nsph - 1; j++)
        {
            if(mu_x[j] < xmax[i] + 5.0 * sigma_x[j] && 
               mu_x[j] > xmin[i] - 5.0 * sigma_x[j] && 
               mu_y[j] < ymax[i] + 5.0 * sigma_y[j] && 
               mu_y[j] > ymin[i] - 5.0 * sigma_y[j] && 
               mu_z[j] < zmax[i] + 5.0 * sigma_z[j] && 
               mu_z[j] > zmin[i] - 5.0 * sigma_z[j])
            {
                total[i] += abs((erf((xmax[i] - mu_x[j]) / sigma_x[j]) - erf((xmin[i] - mu_x[j]) / sigma_x[j])) *
                                (erf((ymax[i] - mu_y[j]) / sigma_y[j]) - erf((ymin[i] - mu_y[j]) / sigma_y[j])) *
                                (erf((zmax[i] - mu_z[j]) / sigma_z[j]) - erf((zmin[i] - mu_z[j]) / sigma_z[j])));
            }
        }
    }

    /* Clean up. */
    Py_XDECREF(xmin_array);
    Py_XDECREF(xmax_array);
    Py_XDECREF(ymin_array);
    Py_XDECREF(ymax_array);
    Py_XDECREF(zmin_array);
    Py_XDECREF(zmax_array);
    Py_XDECREF(mu_x_array);
    Py_XDECREF(mu_y_array);
    Py_XDECREF(mu_z_array);
    Py_XDECREF(sigma_x_array);
    Py_XDECREF(sigma_y_array);
    Py_XDECREF(sigma_z_array);

    return total_array;

}

