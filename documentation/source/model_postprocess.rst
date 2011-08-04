.. _postprocessing:

======================
Post-processing models
======================

You've successfully set up and run your model, but now you want to extract the results and do something useful with them.

Quick look
==========

A convenience script is provided to quickly extract image cubes and physical grids for inspection as FITS files.

::

    # Retrieve all images
    hyperion2fits --images model.rtout

    # Retrieve all physical arrays
    hyperion2fits --physics model.rtout

Ds9 version x.x or later is required to be able to navigate FITS cubes with more than 3 dimensions. This is only meant as a quick look tool. To extract properly scaled and sliced information from the output file, see the sections below.

Extracting Information
======================

The first step to extracting information from the models is to create an instance of the ``Model`` class, giving it the name of the output file minus the extension::

    from hyperion.model import Model
    m = Model('model')

SEDs
----

To extract SEDs, use the ``get_sed`` method::

    wav, nufnu = m.get_sed()

A number of arguments can be passed to ``get_sed()``, for example to select
specific Stokes parameters, inclinations, apertures, to scale the SED to a
specific distance, to convert it to certain units, to extract the SED
originating from different components, etc. For full details about the
available arguments, see the :ref:`getsed` documentation.

What the method returns will depend on the options specified. By default, the
I stokes parameter is returned for all inclinations and apertures. Thus,
``nufnu`` is a data cube with three dimensions (inclinations, apertures, and
wavelengths respectively). If an aperture or an inclination is specified, that
dimension is removed from the array. Thus, specifying both inclination and
aperture makes ``nufnu`` a one-dimensional array.

The default units are microns for ``wav`` and ergs/s for nufnu. If distance is
specified, ``nufnu`` is in ergs/cm^2/s.

If uncertainties are requested, then ``get_sed`` returns three values instead
of two, the third being an uncertainty array with the same dimensions and
units as ``nufnu``::

    wav, nufnu, dnufnu = m.get_sed(uncertainties=True)

Images
------

To extract SEDs, use the ``get_sed`` method::

    wav, nufnu = m.get_image()

As for SEDs, a number of arguments can be passed to ``get_image()``. For full
details about the available arguments, see the :ref:`getimage` documentation.

As for SEDs, the output of the function depends on the options specified. The main difference compared to SEDs is that there are two dimensions for the x and y position in the image instead of the aperture dimension.

Physical Arrays
---------------

Physical arrays (such as specific energy, temperature, density, etc.) can be retrieved using the ``get_physical_grid`` method. For full
details about the available arguments, see the :ref:`getphysics` documentation.

Advanced Topics:
================

.. toctree::
   :maxdepth: 1

   howto_animations.rst
