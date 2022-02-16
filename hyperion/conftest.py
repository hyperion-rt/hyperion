import warnings

from .util.nans import NaNWarning


def pytest_configure(config):
    warnings.simplefilter('error', NaNWarning)
