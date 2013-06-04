from __future__ import print_function, division

import numpy as np

from astropy import log as logger
from ..grid.amr_grid import Grid, Level, AMRGrid


def parse_multi_tuple(string):
    string = string.replace(' ', '')
    string = string.replace(')(', '),(')
    return eval(string)


class Star(object):

    def __init__(self, line):
        values = line.split()
        self.m, self.x, self.y, self.z, self.r, self.mdot, self.burnstate = \
            [float(values[i]) for i in [0, 1, 2, 3, 11, 14, 15]]


class OrionGrid(Grid):

    def __init__(self):
        self.imin, self.imax, self.itype = None, None, None
        self.jmin, self.jmax, self.jtype = None, None, None
        self.kmin, self.kmax, self.ktype = None, None, None
        Grid.__init__(self)

    def read_data(self, filename, offset, quantity_indices, verbose=False):

        if verbose:
            logger.info("Reading %s" % filename)

        gridsize = self.nx * self.ny * self.nz

        f = open(filename, 'rb')
        f.seek(offset)

        header = f.readline().strip()

        p1 = header.find('((') + 2
        p2 = header.find(',', p1)
        n_bytes = int(header[p1:p2])

        p3 = header.find('(', p2) + 1
        p4 = header.find('))', p3)
        bytes = [int(x) for x in header[p3:p4].split()]

        p5 = header.find('(', p4) + 1
        p6 = header.find(',', p5)
        n_bytes = int(header[p5:p6])

        p7 = header.find('(', p6) + 1
        p8 = header.find('))', p7)
        bytes = [int(x) for x in header[p7:p8].split()]

        if bytes == range(1, n_bytes + 1):
            endian = '>'
        elif bytes == range(n_bytes, 0, -1):
            endian = '<'
        else:
            raise Exception("Unexpected bytes: %s" % str(bytes))

        n_components = int(header.split()[-1])

        pos = f.tell()

        for quantity in quantity_indices:

            f.seek(pos + quantity_indices[quantity] * n_bytes * gridsize)
            array = np.fromstring(f.read(n_bytes * gridsize),
                                  dtype='%sf%i' % (endian, n_bytes))
            self.quantities[quantity] = array.reshape(self.nz, self.ny, self.nx)


class OrionLevel(Level):

    def __init__(self):
        self.idxlo = None
        self.idxhi = None
        self.periodicity = None
        self.number = None
        Level.__init__(self)


class OrionAMRGrid(AMRGrid):

    def __init__(self, dirname, quantities, verbose=False, max_level=None):

        self.xmin, self.xmax = None, None
        self.ymin, self.ymax = None, None
        self.zmin, self.zmax = None, None

        AMRGrid.__init__(self)

        # Open file
        f = open('%s/Header' % dirname, 'rb')

        # Read version number
        version = f.readline().strip()

        # Read number of components
        n_quantities = int(f.readline().strip())

        # Read in component names
        available_quantities = [f.readline().strip() for i in range(n_quantities)]

        # If a single quantity is requested as a string, make it into a list
        if isinstance(quantities, basestring):
            if quantities == 'all':
                quantities = available_quantities
            else:
                quantities = [quantities]

        # Make list of wanted quantities, and their indices
        quantity_indices = {}
        for quantity in quantities:
            quantity_indices[quantity] = available_quantities.index(quantity)

        # Read in number of dimensions
        ndim = int(f.readline().strip())
        if ndim != 3:
            raise Exception("Number of dimensions is not 3")

        # Read in time
        creation_time = float(f.readline().strip())

        # Read in maximum level of refinement
        n_levels = int(f.readline().strip()) + 1

        # Create list of levels
        self.levels = [OrionLevel() for i in range(n_levels)]

        if max_level is None:
            max_level = n_levels

        # Read in position of box corners
        self.xmin, self.ymin, self.zmin = [float(x) for x in f.readline().strip().split()]
        self.xmax, self.ymax, self.zmax = [float(x) for x in f.readline().strip().split()]

        # Read in refinement ratios
        refinement_ratios = [int(x) for x in f.readline().strip().split()]

        # Read in next line
        line = f.readline().strip()

        # Split into groups of ndim values
        elements = line.replace(' ', '').replace('((', '(').replace('))', ')')[1:-1].split(')(')
        for level in self.levels:
            level.idxlo = [int(x) for x in elements[3 * i].split(',')]
            level.idxhi = [int(x) for x in elements[3 * i + 1].split(',')]
            level.periodicity = [int(x) for x in elements[3 * i + 2].split(',')]

        # Read in number of steps on each level
        levelsteps = [int(x) for x in f.readline().strip().split()]

        # Read in grid spacing on each level
        # gridspacing = np.zeros((self.ndim, self.maxlevel+1))
        # for level in self.levels:
        #     level.gridspacing = [float(x) for x in f.readline().strip().split()]
        for level in self.levels:
            f.readline()

        # Read in coordinate type
        coordtype = int(f.readline().strip())
        if coordtype != 0:
            raise Exception("coordtype should be zero")

        # Skip dummy line
        f.readline()

        # First part done. Now need to read in individual levels and grids

        # Initialize list of levels

        # Loop through levels
        for level in self.levels[:max_level]:

            level_num, ngrids, creation_time = f.readline().strip().split()
            level.number = int(level_num)
            ngrids = int(ngrids)

            # Initialize grids
            level.grids = [OrionGrid() for igrid in range(ngrids)]

            levelsteps = int(f.readline().strip())

            for grid in level.grids:
                grid.xmin, grid.xmax = [float(x) for x in f.readline().split()]
                grid.ymin, grid.ymax = [float(y) for y in f.readline().split()]
                grid.zmin, grid.zmax = [float(z) for z in f.readline().split()]

            n_quantities_check = 0
            nfiles = 0
            nfilecomp = []

            # Read filename header file
            fname = f.readline().strip()

            fh = open("%s/%s_H" % (dirname, fname))

            fh.readline()
            fh.readline()

            # Read the number of components in multigrid files
            ngridcomp = int(fh.readline())

            if ngridcomp != n_quantities:
                raise Exception("Only some of the components included in multigrid file")

            fh.readline()

            # Read the number of boxes
            ngrids_check = int(fh.readline().strip()[1:].split()[0])

            if ngrids_check != ngrids:
                raise Exception("Number of grids in multigrid file does not match known number")

            # Loop through the grids
            for grid in level.grids:
                values = parse_multi_tuple(fh.readline())
                grid.imin, grid.jmin, grid.kmin = values[0]
                grid.imax, grid.jmax, grid.kmax = values[1]
                grid.itype, grid.jtype, grid.ktype = values[2]
                grid.nx = grid.imax - grid.imin + 1
                grid.ny = grid.jmax - grid.jmin + 1
                grid.nz = grid.kmax - grid.kmin + 1

            fh.readline()
            fh.readline()

            for grid in level.grids:
                string = fh.readline().split(':')[1]
                filename = "%s/Level_%i/%s" % (dirname, level.number, string.split()[0].strip())
                offset = int(string.split()[1])
                grid.read_data(filename, offset, quantity_indices, verbose=verbose)

        # Throw away levels that aren't needed
        self.levels = self.levels[:max_level]


def parse_orion(dirname, quantities='density', verbose=False, max_level=None):

    # Read in grid
    amr_grid = OrionAMRGrid(dirname, quantities=quantities, verbose=verbose, max_level=max_level)

    # Read in star particles
    fs = open('%s/StarParticles' % dirname, 'rb')
    fs.readline()
    stars = []
    for line in fs.readlines():
        stars.append(Star(line))

    return amr_grid, stars
