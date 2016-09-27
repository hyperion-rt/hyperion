import warnings

from astropy.tests.helper import pytest

from .util.nans import NaNWarning



def pytest_addoption(parser):
    parser.addoption('--generate-reference', help="generate reference results for bit-level tests", type="string")
    parser.addoption('--enable-bit-level-tests', help="enable bit-level tests", action="store_true")


def pytest_configure(config):
    warnings.simplefilter('error', NaNWarning)