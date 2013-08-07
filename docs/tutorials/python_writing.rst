Writing files in Python
=======================

Pure Python
-----------

The most basic way to write files in Python is to simply open a file with
write access::

    f = open('file.txt', 'wb')

and to then call the ``write`` method to write to the file::

    f.write("Hello World")

Line returns have to be explicitly included using ``\n``::

    f.write("Line 1\n")
    f.write("line 2\n")

And files should be closed with::

    f.close()

The best way to write out variables with this technique is to use
string formatting which is described in more detail
`here <http://docs.python.org/library/stdtypes.html#string-formatting>`_.
The basic command to format variables into a string is::

    format % variables

where ``format`` is a string containing the format statements and variables is
a tuple of the values, for example::

    >>> print "%s %5.2f %10.4e" % ("name", 3.4, 1.e-10)
    name  3.40 1.0000e-10

We can use this when writing out files, so if we have two lists or arrays of
values ``a`` and ``b`` we can do::

    a = [1,2,3,4,5]
    b = [2,6,4,3,2]

    f = open('file.txt', 'wb')
    for i in range(len(a)):
        f.write("%i %5.2f\n" % (a[i], b[i]))
    f.close()

which will produce a file containing::

    1  2.00
    2  6.00
    3  4.00
    4  3.00
    5  2.00

Numpy
-----

`Numpy`_ provides a function called ``savetxt`` that makes it easy to write out
arrays to files. Given two lists or arrays ``a`` and ``b`` as above, one can
simply do::

   import numpy as np
   a = [1,2,3,4,5]
   b = [2,6,4,3,2]
   np.savetxt('file_numpy.txt', zip(a, b), fmt="%i %5.2f")

which produces exactly the same output as above and avoids the for loop.