#!/usr/bin/env python

import numpy
from setuptools import setup, Extension

ext_modules = [Extension("hyperion.util._integrate_core",
                         ['hyperion/util/_integrate_core.c'],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-Wno-error=declaration-after-statement']),
               Extension("hyperion.util._interpolate_core",
                         ['hyperion/util/_interpolate_core.c'],
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=['-Wno-error=declaration-after-statement']),
               Extension("hyperion.importers._discretize_sph",
                         ['hyperion/importers/_discretize_sph.c'],
                         include_dirs=[numpy.get_include()],
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
                         include_dirs=[numpy.get_include()],
                         extra_compile_args = ['-O2', '-Wno-error=declaration-after-statement'],
                         extra_link_args=['-lstdc++'])]

setup(scripts=['scripts/hyperion', 'scripts/hyperion2fits'],
      ext_modules=ext_modules)
