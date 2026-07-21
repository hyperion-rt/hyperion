#!/usr/bin/env python

import numpy
import sys
from setuptools import setup, Extension

kwargs = {}
kwargs['py_limited_api'] = True
kwargs['include_dirs'] = [numpy.get_include()]
if sys.platform != "win32":
    # GCC 14 (used by manylinux_2_28) promotes these to errors by default; the
    # numpy C-API usage here is fine in practice, so keep them as warnings.
    kwargs['extra_compile_args'] = ['-Wno-error=declaration-after-statement',
                                    '-Wno-error=incompatible-pointer-types']

ext_modules = [Extension("hyperion.util._integrate_core",
                         ['hyperion/util/_integrate_core.c'],
                         **kwargs),
              Extension("hyperion.util._interpolate_core",
                         ['hyperion/util/_interpolate_core.c'],
                         **kwargs),
              Extension("hyperion.importers._discretize_sph",
                         ['hyperion/importers/_discretize_sph.c'],
                         **kwargs),
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
                         extra_link_args=['-lstdc++'],
                         **kwargs)]

setup(scripts=['scripts/hyperion', 'scripts/hyperion2fits'],
      ext_modules=ext_modules)
