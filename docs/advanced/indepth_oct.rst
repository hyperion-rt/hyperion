.. _indepth_oct:

============
Octree Grids
============

`Octrees <http://en.wikipedia.org/wiki/Octree>`_ are hierarchical in nature, and therefore it is easiest to think about setting them up in a recursive manner. To set up an octree, we want to populate a list of booleans (referred to as ``refined``).

The first value indicates whether the parent cell is sub-divided. If it is, then the the second element indicates whether the first cell of the parent cell is sub-divided. If it isn't, then the next value indicates whether the second cell of the parent cell is sub-divided. If it is, then we need to specify the booleans for all the children of that cell before we move to the third cell of the parent cell.

For example, the simplest grid is a single cell that is not sub-divided::

    refined = [False]

The next simplest grid is a single grid cell that is only sub-divided once::

    refined = [True, False, False, False, False, False, False, False, False]

It is easier to picture this as a hierarchy::

    refined = [True,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 ]

If we sub-divide the third sub-cell in the parent cell into cells that are themselves not sub-divided, we get::

    refined = [True,
                 False,
                 False,
                 True,
                   False,
                   False,
                   False,
                   False,
                   False,
                   False,
                   False,
                   False,
                 False,
                 False,
                 False,
                 False,
                 False,
                 ]

and so on. The order of the sub-cells is first along x, then along y, then along z.

In practice, we can make use of a recursive function to set up such a grid. The following example demonstrates this::

    # Set the random seed to make this example predictable
    import random
    random.seed('octree-demo')

    def construct_octree(refined=[True]):

        # Loop over subcells
        for subcell in range(8):

            # Insert criterion for whether cell should be sub-divided. Here we
            # just use a random number to demonstrate.
            divide = random.random() < 0.12

            # Append boolean to overall list
            refined.append(divide)

            # If the cell is sub-divided, recursively divide it further
            if divide:
                construct_octree(refined)

        return refined

    oct = construct_octree()

which gives a refined grid with 65 cells and sub-cells. The length of the list should always be one plus a multiple of 8.


