#!/usr/bin/env python

import os
import sys

if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
builtins._HYPERION_SETUP_ = True

from distutils.core import setup, Extension
from distutils.command.sdist import sdist

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py

from numpy import get_include as get_numpy_include

from hyperion.testing.helper import HyperionTest
from hyperion.version import __version__, __dev__

class custom_sdist(sdist):

    user_options = sdist.user_options + [('unstable', None, "make an unstable release (keep __dev__=True)")]

    def __init__(self, *args, **kwargs):
        sdist.__init__(self, *args, **kwargs)
        self.unstable = False

    def run(self):
        if not self.unstable:
            version_file = 'hyperion/version.py'
            content = open(version_file, 'r').read()
            open(version_file, 'w').write(content.replace('__dev__ = True', "__dev__ = False"))
        try:
            sdist.run(self)
        finally:
            if not self.unstable:
                open(version_file, 'w').write(content)

numpy_includes = get_numpy_include()

cmdclass = {}
cmdclass['build_py'] = build_py
cmdclass['test'] = HyperionTest
cmdclass['sdist'] = custom_sdist

ext_modules = [Extension("hyperion.util._integrate_core",
                         ['hyperion/util/_integrate_core.c'],
                         include_dirs=[numpy_includes]),
               Extension("hyperion.util._interpolate_core",
                         ['hyperion/util/_interpolate_core.c'],
                         include_dirs=[numpy_includes])]

scripts = ['hyperion', 'hyperion2fits']

if __dev__:
    scripts.append('mctherm2hyperion')

setup(name='hyperion',
      version=__version__,
      url='http://www.hyperion-rt.org',
      author='Thomas Robitaille',
      author_email='thomas.robitaille@gmail.com',
      packages=['hyperion',
                'hyperion.conf',
                'hyperion.conf.tests',
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
      package_data={'hyperion.model.tests':['data/*.rtout', 'data/*.hdf5'],
                    'hyperion.testing':['coveragerc']},
      scripts=['scripts/' + x for x in scripts],
      cmdclass=cmdclass,
      ext_modules = ext_modules
     )
