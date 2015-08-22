Convolution with transmission curves
====================================

.. warning:: This feature is currently experimental - use at your own risk!

A new experimental feature has been added to Hyperion in recent releases, which
is the ability to do on-the-fly convolution with spectral transmission curves.
Until this feature was released, the only way to convolve the output from
Hyperion with spectral transmission curves was to output a spectral cube with
high spectral resolution and do the convolution outside of Hyperion. In most
cases this is sufficient, but in cases where the images are high resolution
and/or for multiple viewing angles, the memory requirements increase fast.

It is now possible to instead request that the spectral transmission
convolution is done inside of Hyperion. In Hyperion, we use the term *filter*
to refer to spectral transmission curves, but these can also be total system
spectral transmission curves. Care must be taken when using this feature
because it requires a good understanding of how the input spectral transmission
curves are defined in order to obtain accurate results.

When setting up images using e.g.::

    i = m.add_peeled_images(sed=True, image=False)
    
filters can be added by using::

    f = i.add_filter()
    
instead of doing::

    i.set_wavelength_range(...)
    
It is important to note that the two are incompatible - a given image group can
*either* use filters, *or* a fixed wavelength range. However, multiple filters
can be added to a single image group. The filter properties are then set using
e.g.::

    from astropy import units as u

    f = i.add_filter()
    f.name = 'F2'
    f.spectral_coord = [2, 2.1, 2.2, 2.3, 2.4] * u.micron
    f.transmission = [0., 50, 100, 60, 0.] * u.percent
    f.detector_type = 'energy'
    f.alpha = 1.
    f.central_spectral_coord = 2.15 * u.micron

The attributes to set are the following:

* the ``name`` should simple be a string that can be used later to refer to the
  output in that particular filter.

* the ``spectral_coord`` attribute should be set to the x-axis of the spectral
  transmission curve, and can be in frequency, wavelength, or energy.

* the ``transmission`` attribute should be given as an array or a list of
  values. The absolute values do not matter, because they are re-normalized
  before being used in Hyperion, but the values should give the relative
  transmission as a function of frequency. This should not already be
  multiplied or divided by the frequency. It should simply give at a given
  frequency or wavelength, the relative probability that an energy packet
  will pass through the filter/system.
  
* the ``detector_type`` attribute should be set to either ``'energy'`` or
  ``'photons'``. This is important because for a given input spectral shape for
  the emission, if a detector simply measures photons, then proportionally more
  photons will be detected at longer wavelengths relative to shorter
  wavelengths compared to the ratio of the energy detected at longer
  wavelengths to shorter wavelengths. When using the measurement made with the
  detector, it is therefore important to know whether to take into account this
  bias.

* the ``alpha`` attribute is also related to a subtle issue, which is that when
  measuring a flux through a given filter, what is measured is a total amount
  of energy or photons, but in order to convert this to a monochromatic flux
  :math:`F_\nu`, assumptions need to be made about the underlying spectral
  shape. Examples of this are given in the appendix of `Robitaille et al.
  (2007) <http://adsabs.harvard.edu/abs/2007ApJS..169..328R>`_. The parameter
  ``alpha`` is used to indicate that the underlying spectral shape is
  :math:`\nu^\alpha F_\nu \propto {\rm const}`.

* the ``central_spectral_coord`` attribute gives the spectral coordinate at
  which the monochromatic flux should be given.

Once the filters have been set up, Hyperion runs as usual. The output SEDs and
images will be defined at the central wavelengths of the filters.

Note that filter convolution cannot be used in conjunction with raytracing, nor
with monochromatic radiative transfer.
