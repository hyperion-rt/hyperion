How to efficiently compute pure scattering models
=================================================

In some cases, for example if the wavelength is short enough, it is possible to
ignore dust emission when computing images. In such cases, we can make Hyperion
run faster by disabling the temperature calculation as described in
:ref:`pure-scattering`, and also by producing images only at select wavelengths
using :ref:`monochromatic-rt`. The following model demonstrates how to make an
image of the central region of a flared disk in order to image the inner rim::

    from hyperion.model import AnalyticalYSOModel
    from hyperion.util.constants import au, lsun, rsun, tsun, msun

    # Initialize model
    m = AnalyticalYSOModel()

    # Set up star
    m.star.radius = 1.5 * rsun
    m.star.temperature = tsun
    m.star.luminosity = lsun

    # Set up disk
    d = m.add_flared_disk()
    d.rmin = 10 * rsun
    d.rmax = 30. * au
    d.mass = 0.01 * msun
    d.p = -1
    d.beta = 1.25
    d.r_0 = 10. * au
    d.h_0 = 0.4 * au
    d.dust = 'kmh_lite.hdf5'

    # Set up grid
    m.set_spherical_polar_grid_auto(400, 100, 1)

    # Don't compute temperatures
    m.set_n_initial_iterations(0)

    # Don't re-emit photons
    m.set_kill_on_absorb(True)

    # Use raytracing (only important for source here, since no dust emission)
    m.set_raytracing(True)

    # Compute images using monochromatic radiative transfer
    m.set_monochromatic(True, wavelengths=[1.])

    # Set up image
    i = m.add_peeled_images()
    i.set_image_limits(-13 * rsun, 13 * rsun, -13. * rsun, 13 * rsun)
    i.set_image_size(256, 256)
    i.set_viewing_angles([60.], [20.])
    i.set_wavelength_range(1, 1, 1)

    # Set number of photons
    m.set_n_photons(imaging_sources=10000000, imaging_dust=0,
                    raytracing_sources=100000, raytracing_dust=0)

    # Write out the model and run it in parallel
    m.write('disk.rtin')
    m.run('disk.rtout', mpi=True, n_processes=8)

The dust file is available :download:`here <kmh_lite.hdf5>`. Once this model
has run, we can make a plot of the image (including a linear polarization map)::

    import matplotlib.pyplot as plt
    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc

    mo = ModelOutput('disk.rtout')

    wav, fnu = mo.get_image(inclination=0, units='MJy/sr', distance=300. * pc)
    wav, pol = mo.get_image(inclination=0, stokes='linpol')

    fig = plt.figure(figsize=(8,8))

    # Make total intensity sub-plot

    ax = fig.add_axes([0.1, 0.3, 0.4, 0.4])
    ax.imshow(fnu[:,:,0], extent=[-13, 13, -13, 13],
              interpolation='none', cmap=plt.cm.gist_heat,
              origin='lower', vmin=0., vmax=4e9)
    ax.set_xlim(-13., 13.)
    ax.set_ylim(-13., 13.)
    ax.set_xlabel("x (solar radii)")
    ax.set_ylabel("y (solar radii)")
    ax.set_title("Surface brightness")

    # Make linear polarization sub-plot

    ax = fig.add_axes([0.51, 0.3, 0.4, 0.4])
    im = ax.imshow(pol[:,:,0] * 100., extent=[-13, 13, -13, 13],
                   interpolation='none', cmap=plt.cm.gist_heat,
                   origin='lower', vmin=0., vmax=100.)
    ax.set_xlim(-13., 13.)
    ax.set_ylim(-13., 13.)
    ax.set_xlabel("x (solar radii)")
    ax.set_title("Linear Polarization")
    ax.set_yticklabels('')

    axcb = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    plt.colorbar(im, label="%", cax=axcb)
    fig.savefig('inner_disk.png', bbox_inches='tight')
    
which gives:

.. image:: images/inner_disk.png
   :width: 600px
   :align: center
   
This model takes under 4 minutes to run on 8 cores, which is less than would
normally be required to produce an image with this signal-to-noise in scattered
light.
   

