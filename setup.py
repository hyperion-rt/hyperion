#!/usr/bin/env python

from distutils.core import setup, Extension

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py

from hyperion.testing.helper import HyperionTest

from Cython.Distutils import build_ext

from numpy import get_include as get_numpy_include
numpy_includes = get_numpy_include()

scripts = ['hyperion', 'hyperion2fits', 'mctherm2hyperion']

setup(name='hyperion',
      version='0.8.6',
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
                'hyperion.testing',
                'hyperion.util',
                'hyperion.util.tests'],
      scripts=['scripts/' + x for x in scripts],
      cmdclass={'build_py': build_py, 'test':HyperionTest, 'build_ext':build_ext},
      ext_modules = [Extension("hyperion.util.integrate_core", ['hyperion/util/integrate_core.pyx'], include_dirs=[numpy_includes])],
     )
