=================================
hyperion.grid.CartesianGrid
=================================

.. currentmodule:: hyperion.grid

Base grid class
---------------

.. autoclass:: CartesianGrid

   .. rubric:: Attributes

   .. autosummary::

      ~CartesianGrid.x_wall
      ~CartesianGrid.y_wall
      ~CartesianGrid.z_wall
      ~CartesianGrid.x
      ~CartesianGrid.y
      ~CartesianGrid.z
      ~CartesianGrid.gx
      ~CartesianGrid.gy
      ~CartesianGrid.gz

   .. rubric:: Methods

   .. autosummary::

      ~CartesianGrid.set_walls
      ~CartesianGrid.read
      ~CartesianGrid.read_geometry
      ~CartesianGrid.read_quantities
      ~CartesianGrid.write
      ~CartesianGrid.write_single_array
      ~CartesianGrid.add_derived_quantity

   .. rubric:: Methods (detail)

   .. automethod:: CartesianGrid.set_walls
   .. automethod:: CartesianGrid.read
   .. automethod:: CartesianGrid.read_geometry
   .. automethod:: CartesianGrid.read_quantities
   .. automethod:: CartesianGrid.write
   .. automethod:: CartesianGrid.write_single_array
   .. automethod:: CartesianGrid.add_derived_quantity

.. autoclass:: CartesianGridView

   .. rubric:: Methods

   .. autosummary::

      ~CartesianGridView.append
      ~CartesianGridView.add

   .. rubric:: Methods (detail)

   .. automethod:: CartesianGridView.append
   .. automethod:: CartesianGridView.add

