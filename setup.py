#!/usr/bin/env python

from distutils.core import setup

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py

from hyperion.testing.helper import HyperionTest

scripts = ['ttsre2rtin', 'ttsre2rtin_auto', 'hyperion', 'hyperion2fits', 'mctherm2hyperion']

setup(name='hyperion',
      version='0.8.4',
      packages=['hyperion', 'hyperion.model', 'hyperion.model.tests',
                'hyperion.conf', 'hyperion.densities',
                'hyperion.densities.tests', 'hyperion.dust', 'hyperion.util',
                'hyperion.util.tests', 'hyperion.atmos', 'hyperion.grid',
                'hyperion.grid.tests', 'hyperion.sources',
                'hyperion.sources.tests', 'hyperion.importers',
                'hyperion.testing'],
      scripts=['scripts/' + x for x in scripts],
      cmdclass={'build_py': build_py, 'test':HyperionTest}
     )
