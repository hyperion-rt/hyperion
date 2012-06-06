from __future__ import print_function, division

import os
import tempfile

import h5py
import numpy as np

from ...dust import IsotropicDust
from .. import Model
from ...util.functions import random_id


def random_filename():
    return os.path.join(tempfile.mkdtemp(), random_id())


def get_test_dust():
    dust = IsotropicDust([3.e9, 3.e16], [0.5, 0.5], [1., 1.])
    dust.emissivities.set_lte(dust.optical_properties, n_temp=10, temp_min=0.1, temp_max=1600.)
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
    EXCLUDE_ATTR = ['date_started', 'date_ended', 'cpu_time']

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
    for d in data1:
        a1 = np.array(f1[d])
        a2 = np.array(f2[d])
        assert np.all(a1 == a2)

    # Check if list of groups is the same
    assert group1 == group2

    # Loop over all groups and datasets to check attributes
    for item in ['/'] + data1 + group1:

        # Find all attributes
        attr1 = f1[item].attrs.keys()
        attr1.sort()
        attr2 = f2[item].attrs.keys()
        attr2.sort()
        assert attr1 == attr2

        for a in attr1:
            if a not in EXCLUDE_ATTR:
                assert f1[item].attrs[a] == f2[item].attrs[a]
