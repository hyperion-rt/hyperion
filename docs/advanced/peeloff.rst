Advanced settings for peeled images
===================================

Peeloff origin
--------------

First, it is possible to change the origin relative to which the
peeling-off, and therefore the image extent, are defined. This is set to
``(0., 0., 0.)`` by default, but you can change it using::

    image.set_peeloff_origin((x, y, z))

where ``x``, ``y``, and ``z`` are floating-point values giving the cartesian
coordinates of the peeloff origin. This can be used for example when doing
radiative transfer calculations on a simulation, in order to create images
centered on different sources.

Inside observers
----------------

It is also possible to specify that the images should be calculated from the
point of view of an observer that is inside the coordinate grid, rather than
at infinity. For example, if you are calculating a model of the Milky-Way
Galaxy, you could place yourself at the position of the Sun using::

    from hyperion.util.constants import kpc
    image.set_inside_observer((8.5 * kpc, 0., 0.))

Note that in this case, the peeloff origin cannot be set, and the image
limits, rather than being in cm, should be given in degrees on the sky. Note
also that, like sky coordinates, the x range of the image limits should be
inverted. For example::

    image.set_image_limits(65., -65., -1., 1.)

will produce an image going from l=65 to l=-65, and from b=-1 to b=1. Note
that only images (not SEDs) can be computed for inside observers.

Image depth
-----------

Finally, it is possible to restrict the depth along the line of sight from
which photons should be included, with the default being::

    image.set_depth(-np.inf, np.inf)

The minimum and maximum depth are measured in cm. For standard images, the
depth is taken relative to the plane passing through the origin of the
peeling-off. For images calculated for an observer inside the grid, the
default is::

    image.set_depth(0., np.inf)

where the depth is measured relative and away from the observer.

Note that in this mode, the optical depth used to calculate the peeling off
is the total optical depth to the observer, not just the optical depth in
the slab considered. The slab is only used to determine which emission or
scattering events should be included in the image.
