#!/usr/bin/env python

from distutils.core import setup

scripts = ['ttsre2rtin', 'ttsre2rtin_auto', 'hyperion']

setup(name='hyperion',
      version='0.6.4',
      packages=['hyperion', 'hyperion.model', 'hyperion.conf', 'hyperion.densities',
                'hyperion.dust', 'hyperion.util', 'hyperion.atmos',
                'hyperion.grid', 'hyperion.sources', 'hyperion.importers'],
      package_data={'hyperion': ['data/*.hdf5']},
      scripts=['scripts/' + x for x in scripts])
