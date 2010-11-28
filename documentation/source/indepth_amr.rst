.. _amr_indepth:

=========
AMR Grids
=========

AMR grids are specified by nested objects in the Python, with a layout described in :ref:`model`. These objects can be built in several ways:

Programmatically
================

The following example demonstrates how an AMR grid with 5 fabs and 10 levels per fab can be built programmatically from scratch::

    amr = object()
    amr.levels = []
    for ilevel in range(5):
        level = object()
        level.fabs = []
        for ifab in range(10):
            fab = object()
            fab.xmin, fab.xmax = ..., ...
            fab.ymin, fab.ymax = ..., ...
            fab.zmin, fab.zmax = ..., ...
            fab.nx, fab.ny, fab.nz = 32, 32, 32
            fab.data = np.array(...) # should have shape (nx, ny, nz)
            level.fabs.append(fab)
        amr.levels.append(level)
    
From simulation output
======================

Importing functions are available in ``hyperion.importers`` to convert simulation output to the AMR structure required. At this time, only output from the Orion code can be read in. If the output is contained in a directory  ``directory``, then the AMR structure can be retrieved with::

    from hyperion.importers import parse_orion
    amr = parse_orion('directory')