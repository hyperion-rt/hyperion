Filtering logging messages in Python scripts
============================================

All messages printed by the Hyperion Python routines use the built-in `logging
<https://docs.python.org/3/library/logging.html>`_ module. Each Hyperion module
uses its own logger, named after the module (for example
``hyperion.model.model``), and all of these are children of the ``hyperion``
logger. This means that it is possible to filter messages based on importance.
Messages can have one of several levels:

* ``DEBUG`` (``10``): detailed information, typically of interest only when
  diagnosing problems.

* ``INFO`` (``20``): confirmation that things are working as expected

* ``WARNING`` (``30``): An indication that something unexpected happened, or
  indicative of some problem in the near future. The program is still working
  as expected.

* ``ERROR`` (``40``): due to a more serious problem, the program has not been
  able to perform some function (but no exception is being raised).

* ``CRITICAL`` (``50``): a serious error, indicating that the program itself
  may be unable to continue running (but no exception is being raised).

Note that the ``CRITICAL`` level is unlikely to be used, since critical errors
should raise Exceptions in practice.

Messages which are less severe than the logging threshold are ignored.
Hyperion does not configure logging itself, so the Python defaults apply: only
``WARNING``, ``ERROR`` and ``CRITICAL`` messages are shown, and ``INFO`` and
``DEBUG`` messages are not.

To also show ``INFO`` messages (or ``DEBUG`` messages) in a script, configure
logging at the start of the script, for example::

    import logging
    logging.basicConfig(level=logging.INFO)

This affects messages from all packages that use the ``logging`` module. To
change the threshold for Hyperion messages only, set the level on the
``hyperion`` logger, for example to show all messages including ``DEBUG``::

    import logging
    logging.basicConfig()
    logging.getLogger('hyperion').setLevel(logging.DEBUG)

or to show only ``ERROR`` and ``CRITICAL`` messages from Hyperion::

    import logging
    logging.getLogger('hyperion').setLevel(logging.ERROR)

The level can also be set for a specific part of Hyperion, for example
``logging.getLogger('hyperion.dust')``.

Note that this only applies to messages from the Python routines. The output
of the Fortran code when running a model is not affected by these settings.
