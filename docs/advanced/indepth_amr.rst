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
  the parent level. This applies to both the lower and upper boundaries in
  each direction, so each grid has to cover whole cells of the parent level.
  For example, to refine a single parent cell into eight cells, use one grid
  covering the whole parent cell with ``nx = ny = nz = 2``, rather than eight
  grids with one cell each covering an octant of the parent cell.

If these conditions are not met, then the Fortran Hyperion code will raise an
error.

From simulation output
======================

.. _amr_from_yt:

Using yt
--------

If your simulation can be read in by the `yt <https://yt-project.org/>`_
package (which supports many simulation codes), the easiest way to convert it
to a Hyperion AMR grid is to load it with yt and use the
:meth:`~hyperion.grid.AMRGrid.from_yt` method (this requires yt to be
installed, see :doc:`../installation/installation`)::

    import yt
    from hyperion.grid import AMRGrid

    ds = yt.load('DD0010/moving7_0010')
    amr = AMRGrid.from_yt(ds, quantity_mapping={'density': ('gas', 'density')})

The ``quantity_mapping`` argument is required, and is a dictionary mapping the
name to give each quantity in Hyperion (the key) to the name of the field to
extract from the yt dataset (the value). If a field does not exist, the error
message lists all the fields available in the dataset. Note the following:

* The quantities are converted to g/cm^3, so only density fields can be
  imported in this way. All positions are converted to cm.

* The grid is re-centered so that the center of the yt domain
  (``ds.domain_center``) is at the origin in Hyperion. Apart from this shift,
  the x, y, and z axes are the same in yt and Hyperion, so you should
  subtract ``ds.domain_center`` from any yt position (for example when
  positioning sources).

* Only AMR (patch-based) datasets can be converted, since the conversion uses
  the yt grid hierarchy. Octree and particle-based datasets cannot be imported
  with this method.

The densities should be consistent with the units of the dust opacities (see
:doc:`../important/important`). For example, if the dust opacities are given
per unit dust mass but the yt dataset only contains gas densities, you can
define a dust density field in yt and import that instead::

    def _dust_density(field, data):
        return data[('gas', 'density')].in_units('g/cm**3') * 0.01

    ds.add_field(('gas', 'dust_density'), function=_dust_density,
                 units='g/cm**3', sampling_type='cell')

    amr = AMRGrid.from_yt(ds, quantity_mapping={'density': ('gas', 'dust_density')})

The resulting object can then be used to set up the grid and the density in a
model (see :doc:`../setup/setup_grid`)::

    from hyperion.model import Model

    m = Model()
    m.set_amr_grid(amr)
    m.add_density_grid(amr['density'], 'kmh.hdf5')

Orion
-----

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
    
