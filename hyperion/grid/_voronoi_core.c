#include <Python.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

// Declaration of the voro++ wrapping function.
const char *hyperion_voropp_wrap(int **neighbours, int *max_nn, double **volumes, double **bb_min, double **bb_max,
                  double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                  double const *points, int npoints);

/* Define docstrings */
static char module_docstring[] = "C implementation of utility functions used in Voronoi grids";
static char simplex3d_volume_docstring[] = "Compute the volume of a 3d simplex";
static char neighbours_list_loop_docstring[] = "Inner loop for the computation of the neighbours list";
static char region_in_domain_docstring[] = "Test if region is in domain";
static char voropp_wrapper_docstring[] = "voro++ wrapper";

/* Declare the C functions here. */
static PyObject *_simplex3d_volume(PyObject *self, PyObject *args);
static PyObject *_neighbours_list_loop(PyObject *self, PyObject *args);
static PyObject *_region_in_domain(PyObject *self, PyObject *args);
static PyObject *_voropp_wrapper(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
    {"_simplex3d_volume", _simplex3d_volume, METH_VARARGS, simplex3d_volume_docstring},
    {"_neighbours_list_loop", _neighbours_list_loop, METH_VARARGS, neighbours_list_loop_docstring},
    {"_region_in_domain", _region_in_domain, METH_VARARGS, region_in_domain_docstring},
    {"_voropp_wrapper", _voropp_wrapper, METH_VARARGS, voropp_wrapper_docstring},
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
        Py_XDECREF(array);
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
    Py_XDECREF(array);

    /* Build the output value. */
    PyObject *ret = Py_BuildValue("d", retval);
    if (ret == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Couldn't build output.");
        return NULL;
    }

    return ret;
}

static PyObject *_neighbours_list_loop(PyObject *self, PyObject *args)
{
    PyObject *simplices_obj, *nl_obj;

    if (!PyArg_ParseTuple(args, "OO", &simplices_obj, &nl_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *s_array = PyArray_FROM_OTF(simplices_obj, NPY_INT32, NPY_IN_ARRAY);

    /* Handle invalid input. */
    if (s_array == NULL || !PyList_Check(nl_obj)) {
        PyErr_SetString(PyExc_TypeError, "Invalid input objects.");
        Py_XDECREF(s_array);
        return NULL;
    }

    /* Get the array dimensions. */
    int n = (int)PyArray_DIM(s_array, 0), m = (int)PyArray_DIM(s_array, 1);

    /* Get pointer to the data. */
    int32_t *data = (int32_t *)PyArray_DATA(s_array);

    // Run the loop.
    int i, j1, j2;
    for (i = 0; i < n; ++i) {
        int32_t *simplex = data + i * m;
        for (j1 = 0; j1 < m; ++j1) {
            int32_t s_idx = simplex[j1];
            for (j2 = 0; j2 < m; ++j2) {
                int32_t n_idx = simplex[j2];
                if (n_idx != s_idx) {
                    // This one is a borrowed reference, no need to mess
                    // with refcount.
                    PyObject *tmp_set = PyList_GetItem(nl_obj,s_idx);
                    // This is a new reference.
                    PyObject *n_obj = PyLong_FromLong(n_idx);
                    PySet_Add(tmp_set,n_obj);
                    Py_XDECREF(n_obj);
                }
            }
        }
    }

    // Final cleanup.
    Py_XDECREF(s_array);

    Py_RETURN_NONE;
}

static PyObject *_region_in_domain(PyObject *self, PyObject *args)
{
    PyObject *region_obj, *vertices_obj, *domain_obj;

    if (!PyArg_ParseTuple(args, "OOO", &region_obj, &vertices_obj, &domain_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *v_array = PyArray_FROM_OTF(vertices_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *d_array = PyArray_FROM_OTF(domain_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* Handle invalid input. */
    if (v_array == NULL || d_array == NULL || !PyList_Check(region_obj)) {
        PyErr_SetString(PyExc_TypeError, "Invalid input objects.");
        Py_XDECREF(v_array);
        Py_XDECREF(d_array);
        return NULL;
    }

    double *v_data = (double*)PyArray_DATA(v_array);
    double *d_data = (double*)PyArray_DATA(d_array);

    int ndim = (int)PyArray_DIM(v_array, 1);

    int r_size = PyList_Size(region_obj), i, j;

    for (i = 0; i < r_size; ++i) {
        int vertex_idx = PyLong_AsLong(PyList_GetItem(region_obj,i));
        double *vertex_ptr = v_data + vertex_idx * ndim;
        for (j = 0; j < ndim; ++j) {
            if (vertex_ptr[j] < d_data[2*j] || vertex_ptr[j] > d_data[2*j + 1]) {
                Py_XDECREF(v_array);
                Py_XDECREF(d_array);
                Py_RETURN_FALSE;
            }
        }
    }

    Py_XDECREF(v_array);
    Py_XDECREF(d_array);
    Py_RETURN_TRUE;
}

static PyObject *_voropp_wrapper(PyObject *self, PyObject *args)
{
    PyObject *sites_obj, *domain_obj;

    if (!PyArg_ParseTuple(args, "OO", &sites_obj, &domain_obj))
        return NULL;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *s_array = PyArray_FROM_OTF(sites_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *d_array = PyArray_FROM_OTF(domain_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* Handle invalid input. */
    if (s_array == NULL || d_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid input objects.");
        Py_XDECREF(s_array);
        Py_XDECREF(d_array);
        return NULL;
    }

    double *s_data = (double*)PyArray_DATA(s_array);
    double *d_data = (double*)PyArray_DATA(d_array);

    int nsites = (int)PyArray_DIM(s_array, 0);

    double *volumes = NULL, *bb_min = NULL, *bb_max = NULL;
    int *neighbours = NULL;
    int max_nn;

    // Call the wrapper.
    const char *status = hyperion_voropp_wrap(&neighbours,&max_nn,&volumes,&bb_min,&bb_max,d_data[0],d_data[1],d_data[2],d_data[3],d_data[4],d_data[5],s_data,nsites);

    if (status != NULL) {
        PyErr_SetString(PyExc_RuntimeError, status);
        Py_XDECREF(s_array);
        Py_XDECREF(d_array);
        return NULL;
    }

    // Unfortunately, it seemse like there is no easy way to just re-use the memory we allocated in the wrapper to create a numpy array
    // without copy. See, e.g., here:
    // http://blog.enthought.com/python/numpy-arrays-with-pre-allocated-memory
    // We will just create new numpy arrays and return them for now.
    npy_intp vol_dims[] = {nsites};
    npy_intp neigh_dims[] = {nsites,max_nn};
    npy_intp bb_dims[] = {nsites,3};

    PyObject *vol_array = PyArray_SimpleNew(1,vol_dims,NPY_DOUBLE);
    PyObject *neigh_array = PyArray_SimpleNew(2,neigh_dims,NPY_INT);
    PyObject *bb_min_array = PyArray_SimpleNew(2,bb_dims,NPY_DOUBLE);
    PyObject *bb_max_array = PyArray_SimpleNew(2,bb_dims,NPY_DOUBLE);
    PyObject *max_int = Py_BuildValue("i",INT_MAX);

    if (vol_array == NULL || neigh_array == NULL || bb_min_array == NULL || bb_max_array == NULL || max_int == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
        free(neighbours);
        free(volumes);
        free(bb_min);
        free(bb_max);
        Py_XDECREF(s_array);
        Py_XDECREF(d_array);
        Py_XDECREF(vol_array);
        Py_XDECREF(neigh_array);
        Py_XDECREF(bb_min_array);
        Py_XDECREF(bb_max_array);
        return NULL;
    }

    // Copy over the data.
    memcpy((double*)PyArray_DATA(vol_array),volumes,sizeof(double) * nsites);
    memcpy((int*)PyArray_DATA(neigh_array),neighbours,sizeof(int) * nsites * max_nn);
    memcpy((double*)PyArray_DATA(bb_min_array),bb_min,sizeof(double) * nsites * 3);
    memcpy((double*)PyArray_DATA(bb_max_array),bb_max,sizeof(double) * nsites * 3);

    PyObject *retval = PyTuple_Pack(5,neigh_array,vol_array,bb_min_array,bb_max_array,max_int);

    // Final cleanup.
    free(neighbours);
    free(volumes);
    free(bb_min);
    free(bb_max);
    Py_XDECREF(s_array);
    Py_XDECREF(d_array);
    // NOTE: these need to be cleaned up as PyTuple_Pack will increment the reference count. See:
    // https://mail.python.org/pipermail/capi-sig/2009-February/000222.html
    Py_XDECREF(neigh_array);
    Py_XDECREF(vol_array);
    Py_XDECREF(bb_min_array);
    Py_XDECREF(bb_max_array);
    Py_XDECREF(max_int);

    return retval;
}
