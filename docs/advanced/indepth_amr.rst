.. _amr_indepth:

=========
AMR Grids
=========

AMR grids are specified by nested objects in the Python, with a layout described in :doc:`../setup/setup_grid`. These objects can be built in several ways:

Programmatically
================

The following example demonstrates how an AMR grid can be built programmatically from scratch::

    from hyperion.grid import AMRGrid
    amr = AMRGrid()
    for ilevel in range(nlevels):
        level = amr.add_level()
        for igrid in range(ngrids):
            grid = level.add_grid()
            grid.xmin, grid.xmax = ..., ...
            grid.ymin, grid.ymax = ..., ...
            grid.zmin, grid.zmax = ..., ...
            grid.nx, grid.ny, grid.nz = ..., ..., ...
            grid.quantities['density'] = np.array(...)

where ``nlevels`` is the number of levels in the AMR grid, and ``ngrids`` is the number of grids each each level. The dimensions of the ``np.array(...)`` on the last line should be ``(nz, ny, nx)``.

Requirements
============

The following geometrical requirements have to be respected:

* In a given level, all grids should have the same x, y, and z resolutions,
  i.e. the resolution in each direction can be different, but the resolution
  along a particular direction has to be the same for all grids in the level.

* In a given level, the edges of all grids have to line up with a common grid
  defined by the widths in each direction. For example, if the cell width is
  ``1.0`` in the x direction, one cannot have a grid with ``xmin=0.0`` and one
  with ``xmin=0.5`` since the cell walls in the x direction for these two
  grids do not line up on a common grid.

* The refinement ratio between two levels (the ratio of widths of cells in a
  direction from one level to the next) should be a whole number. The
  refinement ratio can be different for different directions, and can be
  greater than 2.

* The boundaries of grids in a given level have to line up with cell walls in
  the parent level.

If these conditions are not met, then the Fortran Hyperion code will raise an
error.

From simulation output
======================

Importing functions are available in ``hyperion.importers`` to convert simulation output to the AMR structure required. At this time, only output from the Orion code can be read in. If the output is contained in a directory  ``directory``, then the AMR structure can be retrieved with::

    from hyperion.importers import parse_orion
    amr, stars = parse_orion('directory')

The ``stars`` variable is a list of ``Star`` instances. These ``Star`` instances have several attributes, which include:

* ``x``, ``y``, and ``z`` - the position of the star
* ``m``, ``r`` - the mass and radius of the star
* ``mdot`` - the infall rate onto the star

These can be used for example to set up sources of emission in the model::

    # Set up the stars
    for star in stars:
        source = m.add_point_source()
        source.luminosity = lsun
        source.position = (star.x, star.y, star.z)
        source.temperature = 6000.

The above just creates sources with equal temperatures and luminosities, but these can also be set depending on ``m``, ``r``, and ``mdot``.
    
