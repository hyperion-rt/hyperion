Extracting SEDs and Images
==========================

The first step to extracting SEDs and images from the models is to create an instance of the ``ModelOutput`` class, giving it the name of the output file::

    from hyperion.model import ModelOutput
    m = ModelOutput('simple_model.rtout')

SEDs
----

To extract SEDs, use the :meth:`~hyperion.model.ModelOutput.get_sed` method::

    wav, nufnu = m.get_sed()

A number of arguments can be passed to ``get_sed()``, for example to select
specific Stokes parameters, inclinations, apertures, to scale the SED to a
specific distance, to convert it to certain units, to extract the SED
originating from different components, etc. For full details about the
available arguments, see the :meth:`~hyperion.model.ModelOutput.get_sed` documentation.

What the method returns will depend on the options specified. By default, the
I stokes parameter is returned for all inclinations and apertures. Thus,
``nufnu`` is a data cube with three dimensions (inclinations, apertures, and
wavelengths respectively). If an aperture or an inclination is specified, that
dimension is removed from the array. Thus, specifying both inclination and
aperture makes ``nufnu`` a one-dimensional array.

The default units are microns for ``wav`` and ergs/s for nufnu. If distance is
specified, ``nufnu`` is in ergs/cm^2/s.

If uncertainties are requested, then :meth:`~hyperion.model.ModelOutput.get_sed` returns three values instead
of two, the third being an uncertainty array with the same dimensions and
units as ``nufnu``::

    wav, nufnu, dnufnu = m.get_sed(uncertainties=True)

Images
------

To extract images, use the :meth:`~hyperion.model.ModelOutput.get_image` method::

    wav, nufnu = m.get_image()

As for SEDs, a number of arguments can be passed to :meth:`~hyperion.model.ModelOutput.get_image`. For full
details about the available arguments, see the :meth:`~hyperion.model.ModelOutput.get_image` documentation.

As for SEDs, the output of the function depends on the options specified. The main difference compared to SEDs is that there are two dimensions for the x and y position in the image instead of the aperture dimension.
