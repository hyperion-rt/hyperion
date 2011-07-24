#!/usr/bin/env python

from distutils.core import setup

scripts = ['ttsre2rtin', 'ttsre2rtin_auto', 'hyperion', 'hyperion2fits', 'mctherm2hyperion']

setup(name='hyperion',
      version='0.8.3',
      packages=['hyperion', 'hyperion.model', 'hyperion.conf', 'hyperion.densities',
                'hyperion.dust', 'hyperion.util', 'hyperion.atmos',
                'hyperion.grid', 'hyperion.sources', 'hyperion.importers'],
      scripts=['scripts/' + x for x in scripts])
