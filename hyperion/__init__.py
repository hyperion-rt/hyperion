from __future__ import print_function, division

from .version import __version__

try:
    _HYPERION_SETUP_
except NameError:  # we are not currently running setup.py
    try:
        from .util import integrate
    except ImportError:  # indicates the C extension failed to build
        import os
        import sys
        if os.path.exists('setup.py'):
            print("\nYou appear to be importing Hyperion from the source code "
                  "directory. Try changing to a different directory and try "
                  "importing Hyperion again.\n")
            sys.exit(0)
        else:
            raise

# Set up the test function
_test_runner = None


def _get_test_runner():
    from .testing.helper import TestRunner
    return TestRunner(__path__[0])


def test(package=None, test_path=None, args=None, plugins=None,
         verbose=False, pastebin=None, generate_reference=False,
         bit_level_tests=False, coverage=False):
    '''
    Run Hyperion tests using py.test. A proper set of arguments is
    constructed and passed to `pytest.main`.

    Parameters
    ----------
    package : str, optional
        The name of a specific package to test, e.g. 'model' or
        'densities'. If nothing is specified all default Hyperion tests
        are run.

    test_path : str, optional
        Specify location to test by path. May be a single file or
        directory. Must be specified absolutely or relative to the
        calling directory.

    args : str, optional
        Additional arguments to be passed to `pytest.main` in the `args`
        keyword argument.

    plugins : list, optional
        Plugins to be passed to `pytest.main` in the `plugins` keyword
        argument.

    verbose : bool, optional
        Convenience option to turn on verbose output from py.test. Passing
        True is the same as specifying `-v` in `args`.

    pastebin : {'failed','all',None}, optional
        Convenience option for turning on py.test pastebin output. Set to
        'failed' to upload info for failed tests, or 'all' to upload info
        for all tests.

    generate_reference : str
        Generate reference results for bit-level tests

    bit_level_tests : bool
        Run bit-level tests. These are time-consuming tests that check the
        exact validity of the output, but they are disabled by default.

    coverage : bool, optional
        Generate a test coverage report.  The result will be placed in
        the directory htmlcov.

    See Also
    --------
    pytest.main : py.test function wrapped by `run_tests`.

    '''
    test_runner = _get_test_runner()
    return test_runner.run_tests(
        package=package, test_path=test_path, args=args,
        plugins=plugins, verbose=verbose, pastebin=pastebin,
        generate_reference=generate_reference,
        bit_level_tests=bit_level_tests, coverage=coverage)
