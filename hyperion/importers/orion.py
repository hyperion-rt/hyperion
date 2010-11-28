import numpy as np


def parse_multi_tuple(string):
    string = string.replace(' ', '')
    string = string.replace(')(', '),(')
    return eval(string)


class Star(object):

    def __init__(self, line):
        values = line.split()
        self.m, self.x, self.y, self.z, self.r = [float(values[i]) \
                                                  for i in [0, 1, 2, 3, 11]]


class Grid(object):

    def __init__(self):
        self.xmin, self.xmax = None, None
        self.ymin, self.ymax = None, None
        self.zmin, self.zmax = None, None


class Fab(object):

    def __init__(self):
        self.imin, self.imax, self.itype = None, None, None
        self.jmin, self.jmax, self.jtype = None, None, None
        self.kmin, self.kmax, self.ktype = None, None, None
        self.xmin, self.xmax = None, None
        self.ymin, self.ymax = None, None
        self.zmin, self.zmax = None, None
        self.data = {}

    def read_data(self, filename, offset, quantity_index, verbose=False):

        if verbose:
            print "Reading %s" % filename

        fabsize = self.nx * self.ny * self.nz

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

        f.seek(pos + quantity_index * n_bytes * fabsize)
        array = np.fromstring(f.read()[:n_bytes * fabsize],
                              dtype='%sf%i' % (endian, n_bytes))
        self.data = array.reshape(self.nx, self.ny, self.nz)


class Level(object):

    def __init__(self):

        # The level number
        self.level = 0

        # Minimum and maximum position in index space
        self.imin, self.imax = None, None
        self.jmin, self.jmax = None, None
        self.kmin, self.kmax = None, None

        # The grids
        self.grids = []

        # The fabs
        self.fabs = []


class AMR(object):

    def __init__(self, dirname, quantity, verbose=False):

        # Open file
        f = file('%s/Header' % dirname, 'rb')

        # Read version number
        version = f.readline().strip()

        # Read number of components
        n_quantities = int(f.readline().strip())

        # Read in component names
        quantities = [f.readline().strip() for i in range(n_quantities)]

        # Make list of wanted quantities, and their indices
        quantity_index = quantities.index(quantity)

        # Read in number of dimensions
        ndim = int(f.readline().strip())
        if ndim != 3:
            raise Exception("Number of dimensions is not 3")

        # Read in time
        creation_time = float(f.readline().strip())

        # Read in maximum level of refinement
        max_level = int(f.readline().strip())

        # Create list of levels
        self.levels = [Level() for i in range(max_level + 1)]

        # Read in position of box corners
        self.xmin, self.ymin, self.zmin = [float(x) for x in f.readline().strip().split()]
        self.xmax, self.ymax, self.zmax = [float(x) for x in f.readline().strip().split()]

        # Read in refinement ratios
        self.refinement_ratios = [int(x) for x in f.readline().strip().split()]

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
        for level in self.levels:

            level_num, nfabs, creation_time = f.readline().strip().split()
            level.number = int(level_num)
            nfabs = int(nfabs)

            # Initialize grids
            level.fabs = [Fab() for ifab in range(nfabs)]

            levelsteps = int(f.readline().strip())

            for fab in level.fabs:
                fab.xmin, fab.xmax = [float(x) for x in f.readline().split()]
                fab.ymin, fab.ymax = [float(y) for y in f.readline().split()]
                fab.zmin, fab.zmax = [float(z) for z in f.readline().split()]

            n_quantities_check = 0
            nfiles = 0
            nfilecomp = []

            # Read filename header file
            fname = f.readline().strip()

            fh = open("%s/%s_H" % (dirname, fname))

            fh.readline()
            fh.readline()

            # Read the number of components in multifab files
            nfabcomp = int(fh.readline())

            if nfabcomp != n_quantities:
                raise Exception("Only some of the components included in multifab file")

            fh.readline()

            # Read the number of boxes
            nfabs_check = int(fh.readline().strip()[1:].split()[0])

            if nfabs_check != nfabs:
                raise Exception("Number of fabs in multifab file does not match known number")

            # Loop through the fabs
            for fab in level.fabs:
                values = parse_multi_tuple(fh.readline())
                fab.imin, fab.jmin, fab.kmin = values[0]
                fab.imax, fab.jmax, fab.kmax = values[1]
                fab.itype, fab.jtype, fab.ktype = values[2]
                fab.nx = fab.imax - fab.imin + 1
                fab.ny = fab.jmax - fab.jmin + 1
                fab.nz = fab.kmax - fab.kmin + 1

            fh.readline()
            fh.readline()

            for fab in level.fabs:
                string = fh.readline().split(':')[1]
                filename = "%s/Level_%i/%s" % (dirname, level.number, string.split()[0].strip())
                offset = int(string.split()[1])
                fab.read_data(filename, offset, quantity_index, verbose=verbose)

        # Read in star particles
        fs = open('%s/StarParticles' % dirname, 'rb')
        fs.readline()
        self.stars = []
        for line in fs.readlines():
            self.stars.append(Star(line))


def parse_orion(dirname, quantity='density', verbose=False):
    return AMR(dirname, quantity=quantity, verbose=verbose)
