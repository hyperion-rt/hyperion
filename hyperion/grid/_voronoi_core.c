#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <limits.h>
#include <string.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

// Declaration of the voro++ wrapping function.
const char *hyperion_voropp_wrap(int **sparse_neighbours, int **neigh_pos, int *nn, double **volumes, double **bb_min, double **bb_max, double **vertices, int *max_nv,
                  double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double const *points, int npoints, int with_vertices, const char *wall_str, const double *wall_args_arr,
                  int n_wall_args, int with_sampling, int n_samples, double **sample_points, int **sampling_idx, int *tot_samples, int min_cell_samples, int seed, int verbose);

/* Define docstrings */
static char module_docstring[] = "C implementation of utility functions used in Voronoi grids";
static char voropp_wrapper_docstring[] = "voro++ wrapper";

/* Declare the C functions here. */
static PyObject *_voropp_wrapper(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
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

static PyObject *_voropp_wrapper(PyObject *self, PyObject *args)
{
    PyObject *sites_obj, *domain_obj, *wall_args_obj;
    int with_vertices;
    const char *wall_str = "";
    int verbose;
    int with_sampling, n_samples, min_cell_samples, seed;

    if (!PyArg_ParseTuple(args, "OOiiiiii", &sites_obj, &domain_obj, &with_vertices,
        &with_sampling, &n_samples, &min_cell_samples, &seed, &verbose))
    {
        return NULL;
    }

    double wall_args_arr[7];
    int n_wall_args = 0;

    /* Interpret the input objects as `numpy` arrays. */
    PyObject *s_array = PyArray_FROM_OTF(sites_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *d_array = PyArray_FROM_OTF(domain_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

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

    double *volumes = NULL, *bb_min = NULL, *bb_max = NULL, *vertices = NULL, *sample_points = NULL;
    int *sampling_idx = NULL, *sparse_neighbours = NULL, *neigh_pos = NULL;
    int max_nv, tot_samples, nn;

    // Call the wrapper.
    const char *status = hyperion_voropp_wrap(&sparse_neighbours,&neigh_pos,&nn,&volumes,&bb_min,&bb_max,&vertices,&max_nv,
                                              d_data[0],d_data[1],d_data[2],d_data[3],d_data[4],d_data[5],s_data,nsites,with_vertices,
                                              wall_str,wall_args_arr,n_wall_args,with_sampling,n_samples,&sample_points,&sampling_idx,
                                              &tot_samples,min_cell_samples,seed,verbose
                                             );

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
    npy_intp bb_dims[] = {nsites,3};
    npy_intp vert_dims[] = {nsites,max_nv};
    npy_intp spoints_dims[] = {tot_samples,3};
    npy_intp spoints_idx_dims[] = {nsites+1};
    npy_intp sparse_neigh_dims[] = {nn};
    npy_intp neigh_pos_dims[] = {nsites + 1};

    PyObject *vol_array = PyArray_SimpleNew(1,vol_dims,NPY_DOUBLE);
    PyObject *bb_min_array = PyArray_SimpleNew(2,bb_dims,NPY_DOUBLE);
    PyObject *bb_max_array = PyArray_SimpleNew(2,bb_dims,NPY_DOUBLE);
    // NOTE: Py_BuildValue("") is just a safe way to construct None.
    PyObject *vert_array = with_vertices ? PyArray_SimpleNew(2,vert_dims,NPY_DOUBLE) : Py_BuildValue("");
    PyObject *spoints_array = with_sampling ? PyArray_SimpleNew(2,spoints_dims,NPY_DOUBLE) : Py_BuildValue("");
    PyObject *spoints_idx_array = with_sampling ? PyArray_SimpleNew(1,spoints_idx_dims,NPY_INT) : Py_BuildValue("");
    PyObject *sparse_neigh_array = PyArray_SimpleNew(1,sparse_neigh_dims,NPY_INT);
    PyObject *neigh_pos_array = PyArray_SimpleNew(1,neigh_pos_dims,NPY_INT);

    if (vol_array == NULL || bb_min_array == NULL || bb_max_array == NULL ||
        vert_array == NULL || spoints_array == NULL || spoints_idx_array == NULL || sparse_neigh_array == NULL || neigh_pos_array == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "Memory allocation error.");
        free(volumes);
        free(bb_min);
        free(bb_max);
        free(vertices);
        free(sample_points);
        free(sampling_idx);
        free(sparse_neighbours);
        free(neigh_pos);
        Py_XDECREF(s_array);
        Py_XDECREF(d_array);
        Py_XDECREF(vol_array);
        Py_XDECREF(bb_min_array);
        Py_XDECREF(bb_max_array);
        Py_XDECREF(vert_array);
        Py_XDECREF(spoints_array);
        Py_XDECREF(spoints_idx_array);
        Py_XDECREF(sparse_neigh_array);
        Py_XDECREF(neigh_pos_array);
        return NULL;
    }

    // Copy over the data.
    memcpy((double*)PyArray_DATA(vol_array),volumes,sizeof(double) * nsites);
    memcpy((double*)PyArray_DATA(bb_min_array),bb_min,sizeof(double) * nsites * 3);
    memcpy((double*)PyArray_DATA(bb_max_array),bb_max,sizeof(double) * nsites * 3);
    if (with_vertices) {
        memcpy((double*)PyArray_DATA(vert_array),vertices,sizeof(double) * nsites * max_nv);
    }
    if (with_sampling) {
        memcpy((double*)PyArray_DATA(spoints_array),sample_points,sizeof(double) * tot_samples * 3);
        memcpy((int*)PyArray_DATA(spoints_idx_array),sampling_idx,sizeof(int) * (nsites + 1));
    }
    memcpy((int*)PyArray_DATA(sparse_neigh_array),sparse_neighbours,sizeof(int) * nn);
    memcpy((int*)PyArray_DATA(neigh_pos_array),neigh_pos,sizeof(int) * (nsites + 1));

    PyObject *retval = PyTuple_Pack(8,sparse_neigh_array,neigh_pos_array,vol_array,bb_min_array,bb_max_array,vert_array,spoints_array,spoints_idx_array);

    // Final cleanup.
    free(volumes);
    free(bb_min);
    free(bb_max);
    free(vertices);
    free(sample_points);
    free(sampling_idx);
    free(sparse_neighbours);
    free(neigh_pos);
    Py_XDECREF(s_array);
    Py_XDECREF(d_array);
    // NOTE: these need to be cleaned up as PyTuple_Pack will increment the reference count. See:
    // https://mail.python.org/pipermail/capi-sig/2009-February/000222.html
    Py_XDECREF(vol_array);
    Py_XDECREF(bb_min_array);
    Py_XDECREF(bb_max_array);
    Py_XDECREF(vert_array);
    Py_XDECREF(spoints_array);
    Py_XDECREF(spoints_idx_array);
    Py_XDECREF(sparse_neigh_array);
    Py_XDECREF(neigh_pos_array);

    return retval;
}
