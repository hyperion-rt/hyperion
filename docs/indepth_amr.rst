.. _amr_indepth:

=========
AMR Grids
=========

AMR grids are specified by nested objects in the Python, with a layout described in :ref:`model`. These objects can be built in several ways:

Programmatically
================

The following example demonstrates how an AMR grid with 5 grids and 10 levels can be built programmatically from scratch::

    amr = object()
    amr.levels = []
    for ilevel in range(5):
        level = object()
        level.grids = []
        for igrid in range(10):
            grid = object()
            grid.xmin, grid.xmax = ..., ...
            grid.ymin, grid.ymax = ..., ...
            grid.zmin, grid.zmax = ..., ...
            grid.nx, grid.ny, grid.nz = 32, 32, 32
            grid.data = np.array(...) # should have shape (nx, ny, nz)
            level.grids.append(grid)
        amr.levels.append(level)

From simulation output
======================

Importing functions are available in ``hyperion.importers`` to convert simulation output to the AMR structure required. At this time, only output from the Orion code can be read in. If the output is contained in a directory  ``directory``, then the AMR structure can be retrieved with::

    from hyperion.importers import parse_orion
    amr = parse_orion('directory')

As well as a ``levels`` attribute, the amr object retrieved in this way contains a ``stars`` attribute, which is a list of ``Star`` instances. These ``Star`` instances have several attributes, which include:

* ``x``, ``y``, and ``z`` - the position of the star
* ``m``, ``r`` - the mass and radius of the star
* ``mdot`` - the infall rate onto the star

These can be used for example to set up sources of emission in the model::

    # Set up the stars
    for star in amr.stars:
        source = m.add_point_source()
        source.luminosity = lsun
        source.position = (star.x, star.y, star.z)
        source.temperature = 6000.

The above just creates sources with equal temperatures and luminosities, but these can also be set depending on ``m``, ``r``, and ``mdot``.
    