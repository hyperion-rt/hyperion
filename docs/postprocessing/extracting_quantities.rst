Extracting physical quantities
==============================

As described in :doc:`../setup/setup_conf`, it is possible to specify which
gridded physical quantities should be output after the radiative transfer.
For example, by default, the values of the specific energy absorbed in each
cell are output, and this can be used to determine the temperature in each
cell.

To access these gridded physical quantities, first access the model output
file with::

    from hyperion.model import ModelOutput
    m = ModelOutput('output_file.rtout')

then make use of the ``get_quantities`` method to extract a grid object::

    grid = m.get_quantities()

By default, this will extract the physical quantities from the last
iteration, but it is also possible to extract quantities from previous
iterations, e.g.::

    grid = m.get_quantities(iteration=1)

The value of ``iteration`` should be zero-based, so ``1`` indicates the
second iteration.

The ``grid`` variable now contains both the geometrical information about
the grid, and the quantities output in the iteration specified. How this
information is accessed depends on the grid type, as described below.

Regular 3-d grids
-----------------

For regular 3-d grids, the position of the center of the cells can be
accessed with the ``x``, ``y``, and ``z`` attributes for cartesian grids,
the ``w``, ``z``, and ``p`` attributes for cylindrical polar grids, and the
``r``, ``t``, and ``p`` attributes for spherical polar grids. These are 1-d
arrays, but it is also possible to retrieve 3-d versions of these arrays by
adding ``g`` in front, e.g. ``gx``, ``gy``, and ``gz``. It is also possible
to retrieve the original wall positions by adding the ``_wall`` suffix to
the 1-d array names, e.g. ``x_wall``, ``y_wall``, and ``z_wall``.

The physical quantities are stored in a dictionary called ``quantities``,
where each element in the dictionary is a quantity, and each quantity is
stored as a list with as many elements as dust types. Each element is a 3-d
Numpy array. Therefore, you can access for example the ``specific_energy``
values for the first dust type with::

    values = g.quantities['specific_energy'][0]

However, it is also possible to access this with a more convenient
notation::

    values = g['specific_energy'][0].array

Although the temperature values are not ever directly present in the output file, if the specific energy values are present, ``get_quantities`` will automatically calculate and make available the ``temperature`` quantity which can be accessed as above with::

    values = g['temperature'][0].array

In :doc:`../tutorials/tutorial_quantities`, we show how to visualize
this information.

AMR grids
---------

When extracting an AMR grid using ``get_quantities()``, an ``AMRGrid``
object is returned. This object contains an attribute named ``levels`` that
is a list of ``Level`` objects. Each ``Level`` object contains a ``grids``
attribute that is a list of ``Grid`` objects, which in turn have attributes
``xmin``, ``xmax``, ``ymin``, ``ymax``, ``zmin``, ``zmax``, ``nx``, ``ny``,
and ``nz`` which give the boundaries and number of cells in each direction
in the grid (this format is described in more detail
in :doc:`../advanced/indepth_amr`).

Since this is not easy to visualize, Hyperion includes an interface to the
`yt <http://yt-project.org/>`_ package for AMR grids. If you extracted the
quantities with::

    amr = m.get_quantities()

you can call the following method to output a ``StreamStaticOutput`` yt
object that can be directly used for plotting in yt::

    pf = amr.to_yt()

where ``pf`` is a yt ``StaticOutput`` object. See
:doc:`../tutorials/tutorial_quantities_yt` for more details and a plotting
tutorial.

Octree grids
------------

When extracting an Octree grid using ``get_quantities()``, an ``OctreeGrid``
object is returned. The format of this object is described in detail in
:doc:`../advanced/indepth_oct`.

As for AMR grids, Hyperion includes an interface to the `yt
<http://yt-project.org/>`_ package for Octree grids. If you extracted the
quantities with::

    oct = m.get_quantities()

you can call the following method to output a ``StreamStaticOutput`` yt
object that can be directly used for plotting in yt::

    pf = oct.to_yt()

where ``pf`` is a yt ``StaticOutput`` object. See
:doc:`../tutorials/tutorial_quantities_yt` for more details and a plotting
tutorial.

