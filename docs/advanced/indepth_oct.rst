.. _indepth_oct:

============
Octree Grids
============

Setting up an octree grid programmatically
==========================================

`Octrees <http://en.wikipedia.org/wiki/Octree>`_ are hierarchical in nature,
and therefore it is easiest to think about setting them up in a recursive
manner. To set up an octree, we want to populate a list of booleans (referred
to as ``refined``).

The first value indicates whether the parent cell is sub-divided. If it is,
then the the second element indicates whether the first cell of the parent cell
is sub-divided. If it isn't, then the next value indicates whether the second
cell of the parent cell is sub-divided. If it is, then we need to specify the
booleans for all the children of that cell before we move to the third cell of
the parent cell.

For example, the simplest grid is a single cell that is not sub-divided::

    refined = [False]

The next simplest grid is a single grid cell that is only sub-divided once::

    refined = [True, False, False, False, False, False, False, False, False]

It is easier to picture this as a hierarchy::

    refined = [True,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 ]

If we sub-divide the third sub-cell in the parent cell into cells that are themselves not sub-divided, we get::

    refined = [True,
                 False,
                 False,
                 True,
                   False,
                   False,
                   False,
                   False,
                   False,
                   False,
                   False,
                   False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 ]

and so on. The order of the sub-cells is first along x, then along y, then along z.

In practice, we can make use of a recursive function to set up such a grid. The following example demonstrates this::

    # Set the random seed to make this example predictable
    import random
    random.seed('octree-demo')

    def construct_octree(refined=[True]):

        # Loop over subcells
        for subcell in range(8):

            # Insert criterion for whether cell should be sub-divided. Here we
            # just use a random number to demonstrate.
            divide = random.random() < 0.12

            # Append boolean to overall list
            refined.append(divide)

            # If the cell is sub-divided, recursively divide it further
            if divide:
                construct_octree(refined)

        return refined

    oct = construct_octree()

which gives a refined grid with 65 cells and sub-cells. The length of the list should always be one plus a multiple of 8.

Constructing an Octree from a set of SPH particles
==================================================

.. note:: This functionality is currently experimental, so use with care!

Hyperion includes a function, ``construct_octree``, which makes it easy to
produce an Octree for a set of SPH particles, including computing the densities
inside each cell. The function is used as follows::

    from hyperion.importers.sph import construct_octree
    octree = construct_octree(x, y, z, dx, dy, dz, px, py, pz, sigma, mass)

The arguments are:

* the center of the octree grid (``x``, ``y``, ``z``, which should be
  floating-point values)
* the half-width of the octree grid (``dx``, ``dy``, ``dz``, which should be
  floating-point values)
* the positions of the particles (``px``, ``py``, ``pz``, which should be 1-d
  Numpy arrays)
* the sigmas of the Gaussian kernel (``sigma``, which should be a 1-d Numpy
  array)
* the masses of the particles (``mass``, which should be a 1-d Numpy array)

Note that only Gaussian kernels are supported at this time. All values should
be given in cgs. The function returns an :class:`~hyperion.grid.OctreeGrid`
object. This can then be used to set the geometry and the density grid::

    m.set_grid(octree)
    m.add_density_grid(octree['density'][0], dust_file)

A criterion can be specified to halt the refinement of cells. By default, cells
are no longer refined if they contain two or fewer particles. This can be
changed by passing a custom function to ``construct_octree`` using the
``stopping_criterion`` argument. The function passed should take ten
arguments, which are ``x``, ``y``, ``z``, ``dx``, ``dy``, ``dz`` for the
current cell, and the positions and sigmas ``px``, ``py``, ``pz``, and
``sigma`` for the particles in the cell. The function should return ``False``
if the cell should be refined further, and ``True`` otherwise. For example, the
default function can be written as::

    def DEFAULT_STOPPING_CRITERION(x, y, z, dx, dy, dz, px, py, pz, sigma):
        return len(px) <= 2

In addition, the ``construct_octree`` function can take a ``n_levels`` argument
to indicate the maximum number of levels of refinement to allow.

Writing and Reading Octree grids to disk
========================================

Computing an Octree can be a computationally expensive operation, so once it
has been computed, you can write it out to disk using::

    import h5py
    f = h5py.File('my_octree_grid.hdf5', 'w')
    octree.write(f)
    f.close()

You can then read it into a separate script (e.g. the script setting up the
actual model) using::

    import h5py
    from hyperion.grid import OctreeGrid
    f = h5py.File('my_octree_grid.hdf5', 'r')
    octree = OctreeGrid()
    octree.read(f)
    f.close()
