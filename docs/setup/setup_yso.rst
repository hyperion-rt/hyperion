Setting up a YSO
================

Source parameters
-----------------

The stellar luminosity and radius should be set via the following attributes::

    m.star.luminositm = 5 * lsun
    m.star.radius = 2 * rsun

and either the temperature or the spectrum of the source can be set, using::

    m.star.temperature = 10000.

or::

    m.star.spectrum = ...

TODO: spots

Flared disks
------------

Flared disks can be added using the ``add_flared_disk`` method, and capturing the reference to the FlaredDisk object to set the parameters further::

    disk = m.add_flared_disk()
    disk.mass = 0.01 * msun             # Disk mass
    disk.rmin = 10 * rsun               # Inner radius
    disk.rmax = 300 * au                # Outer radius
    disk.h_star = 0.01 * rsun           # Disk scaleheigh at r_star
    disk.alpha = 2.25                   # Radial volume density exponent
    disk.beta = 1.25                    # Disk flaring power

The accretion properties of the disk can be specified in two ways. Either the disk accretion rate can be specified:

    disk.mdot = 1e-6 * msun / yr        # Disk accretion rate

or the accretion luminosity from viscous dissipation:

    disk.lacc = 0.01 * lsun

Note that this accretion luminosity only includes the luminosity down to
``disk.rmin``, and does not include the luminosity from the stellar surface.

Envelopes
---------

Power-law spherically symmetric envelope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest kind of envelope is a spherically symmetric envelope with a power-law distribution in densitm. A power-law envelope can be added using the ``add_powerlaw_envelope`` method, and capturing the reference to the PowerLawEnvelope object to set the parameters further::

    envelope = m.add_powerlaw_envelope()
    envelope.mass = 0.1 * msun          # Envelope mass
    envelope.rmin = au                  # Inner radius
    envelope.rmax = 10000 * au          # Outer radius
    envelope.power = -2                 # Radial power

Ulrich rotationally flattened envelope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A more complex envelope density distribution is that of Ulrich (1976), which consists of a rotationally flattened envelope::

    envelope = m.add_ulrich_envelope()
    envelope.mdot = 1e-4 * msun / yr    # Infall rate
    envelope.rmin = 0.1 * au            # Inner radius
    envelope.rc = 100 * au              # Centrifugal radius
    envelope.rmax = 1000 * au           # Outer radius

Bipolar cavities
----------------

Once an envelope has been created, bipolar cavities can be carved out in it by doing::

    cavitm = envelope.add_bipolar_cavities()
    cavitm.exponent = 1.5               # Shape exponent z~w^exp
    cavitm.r_0 = 1.e-20                 # Radius to specify rho_0 and theta_0
    cavitm.theta_0 = 10                 # Opening angle at r_0 (degrees)
    cavitm.rho_0 = 1.e-20               # Density at r_0
    cavitm.rho_exp = 0.                 # Vertical density exponent

Dust
----

The dust file to use for each component should be specified using the ``dust`` attribute for the component, e.g.::

    disk.dust = 'www003.hdf5'
    envelope.dust = 'kmh.hdf5'
    cavitm.dust = 'kmh_hdf5'

Grid
----

The gridding of the density is done automatically, but the user needs to specify a grid size. Either a spherical polar or cylindrical polar grid can be used. To use the spherical polar grid::

    m.set_spherical_polar_grid_auto(n_r, n_theta, n_phi)

and to use the cylindrical polar grid::

    m.set_cylindrical_polar_grid_auto(n_w, n_z, n_phi)