from __future__ import print_function, division

import os
import tempfile

import h5py
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp

from ...dust import IsotropicDust
from .. import Model
from ...util.functions import random_id


def get_test_dust(set_emissivities=True):
    dust = IsotropicDust([3.e9, 3.e16], [0.5, 0.5], [1., 1.])
    if set_emissivities:
        dust.set_lte_emissivities(n_temp=10, temp_min=0.1, temp_max=1600.)
    return dust


def get_realistic_test_dust():

    nu = [3.e7, 1.e10, 2.e11, 2.e12, 2.e13, 2.e14, 2.e15, 2.e16, 2.e17]
    chi = [1.e-11, 2.e-6, 2.e-3, 0.2, 13., 90., 1000., 700., 700.]
    albedo = [0., 0., 0., 0., 0.1, 0.5, 0.4, 0.4, 0.4]

    dust = IsotropicDust(nu, albedo, chi)
    dust.set_lte_emissivities(n_temp=40, temp_min=0.1, temp_max=100000.)

    return dust


def get_highly_reflective_dust():

    nu = [3.e7, 1.e10, 2.e11, 2.e12, 2.e13, 2.e14, 2.e15, 2.e16, 2.e17]
    chi = np.repeat(100., len(nu))
    albedo = np.repeat(0.7, len(nu))

    dust = IsotropicDust(nu, albedo, chi)
    dust.set_lte_emissivities(n_temp=40, temp_min=0.1, temp_max=100000.)

    return dust


def get_test_model_noimaging():

    model = Model()
    model.set_cartesian_grid([-1., 1.], [-1., 1.], [-1., 1.])
    model.set_n_photons(initial=1, imaging=0)
    model.set_n_initial_iterations(1)

    source = model.add_point_source()
    source.luminosity = 1.
    source.temperature = 1000.

    return model


def assert_identical_results(file1, file2):

    # List of attributes to exclude from checking (time-dependent)
    EXCLUDE_ATTR = ['date_started', 'date_ended', 'cpu_time', 'python_version', 'fortran_version', 'd_min', 'd_max']

    # TODO
    # For now, also exclude 'killed' attributes because they have been moved
    # to a different group, but not worth re-generating all the reference
    # models just for this. However, update this next time the reference
    # models are re-generated.
    EXCLUDE_ATTR += ['killed_photons_geo_initial',
                     'killed_photons_int_initial',
                     'killed_photons_geo',
                     'killed_photons_int']

    # Open both files
    f1 = h5py.File(file1, 'r')
    f2 = h5py.File(file2, 'r')

    # List datasets and groups in file 1
    data1 = []
    group1 = []

    def func(name, obj):
        if isinstance(obj, h5py.Dataset):
            data1.append(name)
        elif isinstance(obj, h5py.Group):
            group1.append(name)
    f1.visititems(func)

     # Visit order is not guaranteed, so sort
    data1.sort()
    group1.sort()

    # List datasets and attributes in file 1
    data2 = []
    group2 = []

    def func(name, obj):
        if isinstance(obj, h5py.Dataset):
            data2.append(name)
        elif isinstance(obj, h5py.Group):
            group2.append(name)
    f2.visititems(func)

     # Visit order is not guaranteed, so sort
    data2.sort()
    group2.sort()

    # Check if list of datasets is the same
    assert data1 == data2

    # Loop over datasets to check content
    TOLERANCE = 1000
    for d in data1:
        a1 = np.array(f1[d])
        a2 = np.array(f2[d])
        if a1.dtype.kind == 'V':  # structured array
            for col in a1.dtype.fields:
                assert_array_almost_equal_nulp(a1[col], a2[col], TOLERANCE)
        else:  # normal array
            assert_array_almost_equal_nulp(a1, a2, TOLERANCE)

    # Check if list of groups is the same
    assert group1 == group2

    # Loop over all groups and datasets to check attributes
    for item in ['/'] + data1 + group1:

        # Find all attributes
        attr1 = f1[item].attrs.keys()
        attr1.sort()
        attr2 = f2[item].attrs.keys()
        attr2.sort()

        for e in EXCLUDE_ATTR:
            if e in attr1:
                attr1.remove(e)
            if e in attr2:
                attr2.remove(e)

        assert attr1 == attr2

        for a in attr1:
            if a not in EXCLUDE_ATTR:
                assert f1[item].attrs[a] == f2[item].attrs[a]
