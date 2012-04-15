#!/usr/bin/env python

from distutils.core import setup

scripts = ['ttsre2rtin', 'ttsre2rtin_auto', 'hyperion', 'hyperion2fits', 'mctherm2hyperion']

setup(name='hyperion',
      version='0.8.4',
      packages=['hyperion', 'hyperion.model', 'hyperion.model.tests',
                'hyperion.conf', 'hyperion.densities',
                'hyperion.densities.tests', 'hyperion.dust', 'hyperion.util',
                'hyperion.util.tests', 'hyperion.atmos', 'hyperion.grid',
                'hyperion.grid.tests', 'hyperion.sources',
                'hyperion.sources.tests', 'hyperion.importers'],
      scripts=['scripts/' + x for x in scripts])
