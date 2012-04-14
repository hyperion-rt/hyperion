Filtering logging messages in Python scripts
============================================

All messages printed by the Hyperion Python routines use the built-in `logging
<http://docs.python.org/library/logging.html>`_ module. This means that it is
possible to filter messages based on importance. Messages can have one of
several levels:

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

It is possible to specify a threshold for logging messages. Messages which are
less severe than this threshold will be ignored. The default threshold in
Hyperion is ``20`` (``INFO``), indicating that all the above messages will be
shown except ``DEBUG``. Using 40 for example would cause only ``ERROR`` and
``CRITICAL`` messages to be shown.

By directly accessing the Hyperion logger
-----------------------------------------

If you want to filter different messages in different scripts, you can
directly access the logger and set the level manually::

    from hyperion.util.logger import logger
    logger.setLevel(10)

Using the Hyperion configuration file
-------------------------------------

If you want to always filter the same messages for all projects on a given
computer, you can create a ``.hyperionrc`` file in your home directory,
containing the following entries related to logging::

    [logging]
    color:yes
    level:0

