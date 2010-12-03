.. _nonrootpython:

============================
Non-root Python installation
============================

If you do not have root access to the computer you wish to use Hyperion on, you will need to use the following instructions to allow packages to be installed in your home directory automatically.

First, create a file called ``.pydistutils.cfg`` in your home directory, with the following contents::

    [install]
    install_lib = ~/usr/python/$py_version_short/site-packages
    install_scripts = ~/usr/bin
    install_data = ~/usr/share
    
Then, add::

    export PYTHONPATH=$HOME/usr/python/2.6/site-packages/:$PYTHONPATH
    export PATH=$HOME/usr/bin:$PATH

to your ``.bashrc`` file (or the equivalent ``setenv`` commands if you are using ``csh``).

Then, close the terminal and open a new one for the changes to take effect. Any package you install using the standard ``python setup.py install`` will now be placed in ``~/usr`` and will be picked up by Python.