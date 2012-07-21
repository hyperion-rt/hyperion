#!/usr/bin/env python

import os

from distutils.core import setup, Extension
from distutils.command import sdist

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py

from numpy import get_include as get_numpy_include

from hyperion.testing.helper import HyperionTest
from hyperion.version import __version__, __dev__

numpy_includes = get_numpy_include()

cmdclass = {}
cmdclass['build_py'] = build_py
cmdclass['test'] = HyperionTest
cmdclass['sdist'] = sdist.sdist

if __dev__:
    from Cython.Distutils import build_ext
    ext_modules = [Extension("hyperion.util.integrate_core",
                             ['hyperion/util/integrate_core.pyx'],
                             include_dirs=[numpy_includes]),
                   Extension("hyperion.util.interpolate_core",
                             ['hyperion/util/interpolate_core.pyx'],
                             include_dirs=[numpy_includes])]
    cmdclass['build_ext'] = build_ext
else:
    ext_modules = [Extension("hyperion.util.integrate_core",
                             ['hyperion/util/integrate_core.c'],
                             include_dirs=[numpy_includes]),
                   Extension("hyperion.util.interpolate_core",
                             ['hyperion/util/interpolate_core.c'],
                             include_dirs=[numpy_includes])]

scripts = ['hyperion', 'hyperion2fits', 'mctherm2hyperion']

setup(name='hyperion',
      version=__version__,
      packages=['hyperion',
                'hyperion.conf',
                'hyperion.densities',
                'hyperion.densities.tests',
                'hyperion.dust',
                'hyperion.dust.tests',
                'hyperion.grid',
                'hyperion.grid.tests',
                'hyperion.importers',
                'hyperion.model',
                'hyperion.model.tests',
                'hyperion.sources',
                'hyperion.sources.tests',
                'hyperion.sphinx',
                'hyperion.sphinx.ext',
                'hyperion.testing',
                'hyperion.util',
                'hyperion.util.tests'],
      package_data={'hyperion.model.tests':['data/*.rtout', 'data/*.hdf5']},
      scripts=['scripts/' + x for x in scripts],
      cmdclass=cmdclass,
      ext_modules = ext_modules
     )
