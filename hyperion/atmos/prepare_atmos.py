from __future__ import print_function, division

import os

import atpy
import numpy as np

from ..util.constants import c
from ..util.logger import logger


def prepare_atmos(filename, output=None):

    logger.info("Converting atmosphere file %s" % os.path.basename(filename))

    atmos = np.loadtxt(filename, comments='#', \
                dtype=[('wav', float), ('fnu', float)])

    t = atpy.Table()
    t.add_column('nu', c * 1.e4 / atmos['wav'])
    t.add_column('wav', atmos['wav'])
    t.add_column('fnu', atmos['fnu'])
    t.sort('nu')

    if output:
        print("-> Writing %s" % os.path.basename(output))
        t.write(output, overwrite=True)
    else:
        return t
