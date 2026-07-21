import warnings

# Use a non-interactive matplotlib backend for the test suite so that plotting
# tests do not require a GUI toolkit (e.g. they would otherwise fail on Windows
# CI where Tcl/Tk is not set up). This must happen before pyplot is imported.
import matplotlib
matplotlib.use('Agg')

from .util.nans import NaNWarning


def pytest_configure(config):
    warnings.simplefilter('error', NaNWarning)
