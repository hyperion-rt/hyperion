An introduction to Matplotlib
=============================

Before we start plotting SEDs and images, let's see how Matplotlib works. The
first thing we want to do is to import Matplotlib::

    import matplotlib.pyplot as plt

In general, if you are plotting from a script as opposed to interactively, you
probably don't want each figure to pop up to the screen before being written
out to a file. In that case, use::

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

This sets the default backend for plotting to ``Agg``, which is a
non-interactive backend.

We now want to create a figure, which we do using::

    fig = plt.figure()

This creates a new figure, which we can now access via the ``fig`` variable
(``fig`` is just the variable name we chose, you can use anything you want).
The figure is just the container for anything you are going to plot, so we
next need to add a set of axes. The simplest way to do this is::

    ax = fig.add_subplot(1,1,1)

The arguments for ``add_subplot`` are the number of subplots in the x and y
direction, followed by the index of the current one.

The most basic command you can now type is ``plot``, which draws a line plot.
The basic arguments are the x and y coordinates of the vertices::

    ax.plot([1,2,3,4],[2,3,2,1])

Now we can simply save the plot to a figure. A number of file types are
supported, including PNG, JPG, EPS, and PDF. Let's write our masterpiece out
to a PNG file::

    fig.savefig('line_plot.png')

.. note:: If you are using a Mac, then writing out in PNG while you are
          working on the plot is a good idea, because if you open it in
          Preview.app, it will update automatically whenever you run
          the script (you just need to click on the plot window to make
          it update). When you are happy with your plot, you can always
          switch to EPS.

Our script now looks like::

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot([1,2,3,4],[2,3,2,1])
    fig.savefig('line_plot.png')