Extracting SEDs and Images
==========================

The first step to extracting SEDs and images from the models is to create an instance of the ``ModelOutput`` class, giving it the name of the output file::

    from hyperion.model import ModelOutput
    m = ModelOutput('simple_model.rtout')

SEDs
----

To extract SEDs, use the :meth:`~hyperion.model.ModelOutput.get_sed` method::

    sed = m.get_sed()

A number of arguments can be passed to
:meth:`~hyperion.model.ModelOutput.get_sed`, for example to select specific
Stokes parameters, inclinations, apertures, to scale the SED to a specific
distance, to convert it to certain units, to extract the SED originating from
different components, etc. For full details about the available arguments, see
the :meth:`~hyperion.model.ModelOutput.get_sed` documentation. The method
returns a single :class:`~hyperion.model.SED` object that contains e.g. the
wavelengths (``sed.wav``), frequencies (``sed.nu``), values (i.e. fluxes, flux
densities, or polarization values; ``sed.val``), and optionally uncertainties
(``sed.unc``). See :class:`~hyperion.model.SED` for the full list of the
available attributes.

By default, the I stokes parameter is returned for all inclinations and
apertures, and ``sed.val`` is a data cube with three dimensions (inclinations,
apertures, and wavelengths respectively). If an aperture or an inclination is
specified, that dimension is removed from the array. Thus, specifying both
inclination and aperture makes ``sed.val`` a one-dimensional array.

The default units are microns for ``sed.wav`` and ergs/s for ``sed.val``. If a
distance is specified when extracting the SED, ``sed.val`` is in ergs/cm^2/s
by default.

If uncertainties are requested, then ``sed.unc`` is set, which is uncertainty
array with the same dimensions and units as ``sed.val``::

    sed = m.get_sed(uncertainties=True)

See :doc:`../tutorials/tutorial_seds` for an example of extracting SEDs from a
model.

Images
------

To extract images, use the :meth:`~hyperion.model.ModelOutput.get_image`
method::

    image = m.get_image()

Similarly to SEDs, a number of arguments can be passed to
:meth:`~hyperion.model.ModelOutput.get_image`. For full details about the
available arguments, see the :meth:`~hyperion.model.ModelOutput.get_image`
documentation. This method returns a single :class:`~hyperion.model.Image`
object that contains e.g. the wavelengths (``image.wav``), frequencies
(``image.nu``), values (i.e. fluxes, flux densities, or polarization
values; ``image.val``), and optionally uncertainties (``image.unc``). See
:class:`~hyperion.model.Image` for the full list of the available attributes.

As for SEDs, the attributes of the image will depend on the options specified.
The main difference compared to SEDs is that there are two dimensions for the x
and y position in the image instead of the aperture dimension.

See :doc:`../tutorials/tutorial_images` for an example of extracting images
from a model.

