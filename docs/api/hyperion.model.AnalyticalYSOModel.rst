=================================
hyperion.model.AnalyticalYSOModel
=================================

:class:`~hyperion.model.AnalyticalYSOModel` inherits from :class:`~hyperion.model.Model`, so all methods and attributes present in the latter can be used with the former, but they are not listed here.

.. currentmodule:: hyperion.model

.. autoclass:: AnalyticalYSOModel

   .. rubric:: Adding density structures

   .. autosummary::

      ~AnalyticalYSOModel.add_flared_disk
      ~AnalyticalYSOModel.add_alpha_disk
      ~AnalyticalYSOModel.add_power_law_envelope
      ~AnalyticalYSOModel.add_ulrich_envelope
      ~AnalyticalYSOModel.add_ambient_medium

   .. rubric:: Setting the grid automatically

   .. autosummary::

      ~AnalyticalYSOModel.set_spherical_polar_grid_auto
      ~AnalyticalYSOModel.set_cylindrical_polar_grid_auto

   .. rubric:: Miscellaneous

   .. autosummary::

      ~AnalyticalYSOModel.evaluate_optically_thin_radii
      ~AnalyticalYSOModel.setup_magnetospheric_accretion

   .. rubric:: Methods (detail)

   .. automethod:: AnalyticalYSOModel.add_flared_disk
   .. automethod:: AnalyticalYSOModel.add_alpha_disk
   .. automethod:: AnalyticalYSOModel.add_power_law_envelope
   .. automethod:: AnalyticalYSOModel.add_ulrich_envelope
   .. automethod:: AnalyticalYSOModel.add_ambient_medium

   .. automethod:: AnalyticalYSOModel.set_spherical_polar_grid_auto
   .. automethod:: AnalyticalYSOModel.set_cylindrical_polar_grid_auto

   .. automethod:: AnalyticalYSOModel.evaluate_optically_thin_radii
   .. automethod:: AnalyticalYSOModel.setup_magnetospheric_accretion
