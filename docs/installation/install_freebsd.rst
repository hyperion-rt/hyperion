Installing dependencies on FreeBSD (requires root)
--------------------------------------------------
I created this guide after getting Hyperion installed on FreeBSD 10.2 with Python 3.

Fortran dependencies
^^^^^^^^^^^^^^^^^^^^

The Fortran Hyperion code requires `HDF5 <http://www.hdfgroup.org/HDF5/>`_ and `MPICH2 <http://www.mpich.org/>`_.


MPICH2
::::::

You can install MPICH2 with the package management systen::

    pkg install mpich2
    
Copy the mpd.conf file to the user home directory that you will be running the simulations from and chown it accordingly. (I hope you followed the note of the installer and changed the secret in it!)::
    
    cp /usr/local/etc/mpd.conf /home/<user>/.mpd.conf
    chown <user>:<user> /home/<user>/.mpd.conf

Run the following command as that user::
    
    mpdboot
    
HDF5
::::::

Since the precompiled version of HDF5 that can be installed using `pkg` doesn't contain the Fortran bindings we have to use the ports collection.::

    cd /usr/ports/science/hdf5
    make config

Make sure to select the Fortran language support!::

    make

This will take some time to compile everything necessary. Once done type the following to install it.::

    make install


Python dependencies
^^^^^^^^^^^^^^^^^^^

I used Python3::
    
    pkg install python3
    
Since Pip isn't installed initially I had to run the command::
    
    python3 -m ensurepip
    
Now pip can be used to install the remaining python dependencies by being called with `pip3`.
To make sure that pkg-config and png is installed::
    
    pkg install pkgconf png
    
After this just go on installing the necessary Python dependencies using pip(3).
    
Hyperion
^^^^^^^^

You are now ready to install Hyperion. Follow the guide of how to install and compile Hyperion.