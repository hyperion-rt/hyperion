Computing dust properties
=========================

A standard algorithm used to compute scattering and absorption properties by
spherical dust particles is that provided by `Bohren and Huffman (1983)`_ (also
known as ``bhmie``). As part of the Hyperion development, we have written a
wrapper around the improved version of `bhmie.f`_ developed by Bruce Draine
that makes it easy to compute the dust properties for mixes of dust
compositions and grain sizes. The program can be found
`here <https://github.com/hyperion-rt/bhmie>`_ along with instructions for
installation and usage.

Note that this code is only able to compute properties for simple spherical
grains with no mantles.

 .. admonition:: Disclaimer

    We cannot guarantee that the wrapper (or the original ``bhmie`` code) are
    bug-free - you should sign up to the `mailing list
    <https://groups.google.com/forum/?fromgroups#!forum/hyperion-announce>`_ to
    ensure that you are informed as soon as bugs are identified and fixed.

In order to be able to easily read computed properties into Hyperion, the
easiest way is to set the ``format`` parameter to ``2``. Then, if the
``prefix`` and ``format`` are set to e.g.::

    'my_dust' = prefix for results
    2 = output format (see README)
    
You can create a Hyperion dust object with::

    from hyperion.dust import BHDust
    d = BHDust('my_dust')
    
and as described in :doc:`../setup/setup_dust`, you can visualize and write out
the dust object with::

    d.plot('my_dust.png')
    d.write('my_dust.hdf5')

Please `let us know`_ if you run into any issues!

.. _bhmie: https://github.com/hyperion-rt/bhmie

.. _bhmie.f: http://www.astro.princeton.edu/~draine/scattering.html

.. _Bohren and Huffman (1983): http://adsabs.harvard.edu/abs/1983asls.book.....B

.. _let us know: http://www.github.com/hyperion-rt/hyperion/issues
