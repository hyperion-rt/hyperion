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
                         include_dirs=[numpy_includes],
                         extra_compile_args=['-Wno-error=declaration-after-statement']),
               Extension("hyperion.util._interpolate_core",
                         ['hyperion/util/_interpolate_core.c'],
                         include_dirs=[numpy_includes],
                         extra_compile_args=['-Wno-error=declaration-after-statement']),
               Extension("hyperion.importers._discretize_sph",
                         ['hyperion/importers/_discretize_sph.c'],
                         include_dirs=[numpy_includes],
                         extra_compile_args=['-Wno-error=declaration-after-statement']),
               Extension("hyperion.grid._voronoi_core",
                         ['hyperion/grid/_voronoi_core.c',
                          'hyperion/grid/voropp_wrap.cc',
                          'hyperion/grid/voro++/c_loops.cc',
                          'hyperion/grid/voro++/cell.cc',
                          'hyperion/grid/voro++/common.cc',
                          'hyperion/grid/voro++/container.cc',
                          'hyperion/grid/voro++/container_prd.cc',
                          'hyperion/grid/voro++/pre_container.cc',
                          'hyperion/grid/voro++/unitcell.cc',
                          'hyperion/grid/voro++/v_base.cc',
                          'hyperion/grid/voro++/v_compute.cc',
                          'hyperion/grid/voro++/wall.cc'],
                         extra_compile_args = ['-O2', '-Wno-error=declaration-after-statement'],
                         include_dirs=[numpy_includes])]

scripts = ['hyperion', 'hyperion2fits']

if __dev__:
    scripts.append('mctherm2hyperion')

setup(name='Hyperion',
      version=__version__,
      url='http://www.hyperion-rt.org',
      description='Monte-Carlo Radiative Transfer Code',
      long_description='Monte-Carlo Radiative Transfer Code',
      author='Thomas Robitaille',
      author_email='thomas.robitaille@gmail.com',
      license='BSD',
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
                'hyperion.importers.tests',
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
                    'hyperion.importers.tests':['data/*.hdf5'],
                    'hyperion.grid.tests':['data/*.hdf5'],
                    'hyperion.testing':['coveragerc']},
      scripts=['scripts/' + x for x in scripts],
      cmdclass=cmdclass,
      ext_modules = ext_modules
     )
