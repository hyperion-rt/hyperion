#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include <math.h>

/* Define docstrings */
static char module_docstring[] = "Helpers for discretizing SPH particles";
static char discretize_sph_docstring[] = "Discretize SPH particles onto given cells";
static char get_positions_widths_docstring[] = "Get positions and widths of all cells";

/* Declare the C functions here. */
static PyObject *_discretize_sph_func(PyObject *self, PyObject *args);
static PyObject *_get_positions_widths(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
    {"_discretize_sph_func", _discretize_sph_func, METH_VARARGS, discretize_sph_docstring},
    {"_get_positions_widths", _get_positions_widths, METH_VARARGS, get_positions_widths_docstring},
    {NULL, NULL, 0, NULL}
};

int recursive_position_width(int i, long *refined,
                             double x, double y, double z,
                             double dx, double dy, double dz,
                             double *xc, double *yc, double *zc,
                             double *xw, double *yw, double *zw);

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

    PyObject *mu_x_obj;
    PyObject *mu_y_obj;
    PyObject *mu_z_obj;
    PyObject *sigma_obj;
    PyObject *mass_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOO", &xmin_obj, &xmax_obj,
                                                &ymin_obj, &ymax_obj,
                                                &zmin_obj, &zmax_obj,
                                                &mu_x_obj, &mu_y_obj, &mu_z_obj,
                                                &sigma_obj, &mass_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *xmin_array = PyArray_FROM_OTF(xmin_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *xmax_array = PyArray_FROM_OTF(xmax_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *ymin_array = PyArray_FROM_OTF(ymin_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *ymax_array = PyArray_FROM_OTF(ymax_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *zmin_array = PyArray_FROM_OTF(zmin_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *zmax_array = PyArray_FROM_OTF(zmax_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *mu_x_array = PyArray_FROM_OTF(mu_x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *mu_y_array = PyArray_FROM_OTF(mu_y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *mu_z_array = PyArray_FROM_OTF(mu_z_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *sigma_array = PyArray_FROM_OTF(sigma_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *mass_array = PyArray_FROM_OTF(mass_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (xmin_array == NULL || xmax_array == NULL ||
        ymin_array == NULL || ymax_array == NULL ||
        zmin_array == NULL || zmax_array == NULL ||
        mu_x_array == NULL || mu_y_array == NULL || mu_z_array == NULL ||
        sigma_array == NULL || mass_array == NULL) {
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
        Py_XDECREF(sigma_array);
        Py_XDECREF(mass_array);
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
        nsph != (int)PyArray_DIM(sigma_array, 0) ||
        nsph != (int)PyArray_DIM(mass_array, 0)
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
        Py_XDECREF(sigma_array);
        Py_XDECREF(mass_array);
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
        Py_XDECREF(sigma_array);
        Py_XDECREF(mass_array);
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
    double *sigma = (double*)PyArray_DATA(sigma_array);
    double *mass = (double*)PyArray_DATA(mass_array);
    double *total = (double*)PyArray_DATA(total_array);

    /* Calculate total */

    int i, j;
    double norm;

    /* Loop over all cells */
    for (i = 0; i < ncells; i++)
    {

        total[i] = 0.;

        /* Loop over all SPH particles */

        for (j = 0; j < nsph; j++)
        {
            if(mu_x[j] < xmax[i] + 3.0 * sigma[j] &&
               mu_x[j] > xmin[i] - 3.0 * sigma[j] &&
               mu_y[j] < ymax[i] + 3.0 * sigma[j] &&
               mu_y[j] > ymin[i] - 3.0 * sigma[j] &&
               mu_z[j] < zmax[i] + 3.0 * sigma[j] &&
               mu_z[j] > zmin[i] - 3.0 * sigma[j])
            {
                norm = 1. / (sigma[j] * sqrt(2.));
                total[i] += fabs((erf((xmax[i] - mu_x[j]) * norm) - erf((xmin[i] - mu_x[j]) * norm)) *
                                 (erf((ymax[i] - mu_y[j]) * norm) - erf((ymin[i] - mu_y[j]) * norm)) *
                                 (erf((zmax[i] - mu_z[j]) * norm) - erf((zmin[i] - mu_z[j]) * norm))) * 0.125 * mass[j];
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
    Py_XDECREF(sigma_array);
    Py_XDECREF(mass_array);

    return total_array;

}

static PyObject *_get_positions_widths(PyObject *self, PyObject *args)
{

    PyObject *refined_obj;
    double x0, y0, z0;
    double dx, dy, dz;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "Odddddd", &refined_obj, &x0, &y0, &z0, &dx, &dy, &dz))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *refined_array = PyArray_FROM_OTF(refined_obj, NPY_LONGLONG, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an `Exception`. */
    if (refined_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't parse the input arrays.");
        Py_XDECREF(refined_array);
        return NULL;
    }

    /* How many cells are there? */
    int ncells = (int)PyArray_DIM(refined_array, 0);

    /* Build the output arrays */

    npy_intp dims[1];
    dims[0] = ncells;

    PyObject *xc_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (xc_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_XDECREF(refined_array);
        return NULL;
    }

    PyObject *yc_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (yc_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_XDECREF(refined_array);
        return NULL;
    }

    PyObject *zc_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (zc_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_XDECREF(refined_array);
        return NULL;
    }

    PyObject *xw_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (xw_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_XDECREF(refined_array);
        return NULL;
    }

    PyObject *yw_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (yw_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_XDECREF(refined_array);
        return NULL;
    }

    PyObject *zw_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (zw_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Couldn't build output array");
        Py_XDECREF(refined_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    long *refined = (long*)PyArray_DATA(refined_array);
    double *xc = (double*)PyArray_DATA(xc_array);
    double *yc = (double*)PyArray_DATA(yc_array);
    double *zc = (double*)PyArray_DATA(zc_array);
    double *xw = (double*)PyArray_DATA(xw_array);
    double *yw = (double*)PyArray_DATA(yw_array);
    double *zw = (double*)PyArray_DATA(zw_array);

    /* Compute cell properties */
    int i;
    i = recursive_position_width(0, refined, x0, y0, z0, dx, dy, dz, xc, yc, zc, xw, yw, zw);

    if(i != ncells - 1) {
        PyErr_SetString(PyExc_TypeError, "An error occurred when retrieving the cell properties");
    }

    // return xc_array, yc_array, zc_array;
    return Py_BuildValue("OOOOOO", xc_array, yc_array, zc_array, xw_array, yw_array, zw_array);

}


int recursive_position_width(int i, long *refined,
                             double x, double y, double z,
                             double dx, double dy, double dz,
                             double *xc, double *yc, double *zc,
                             double *xw, double *yw, double *zw) {

    int ix, iy, iz;

    xc[i] = x;
    yc[i] = y;
    zc[i] = z;
    xw[i] = dx;
    yw[i] = dy;
    zw[i] = dz;

    if(refined[i] == 1) {

        dx *= 0.5;
        dy *= 0.5;
        dz *= 0.5;

        for(iz=-1;iz<2;iz+=2) {
            for(iy=-1;iy<2;iy+=2) {
                for(ix=-1;ix<2;ix+=2) {
                    i = recursive_position_width(i + 1, refined,
                                                 x + ix * dx, y + iy * dy, z + iz * dz,
                                                 dx, dy, dz,
                                                 xc, yc, zc, xw, yw, zw);
                }
            }
        }

    }

    return i;

}
