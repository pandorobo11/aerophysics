API and references
==================

The API reference is organized by calculation domain. Public objects are
documented once, in their defining module. Frequently used classes and
functions are also re-exported from ``aerophysics`` for concise imports:

.. code-block:: python

   from aerophysics import FlightCondition, oblique_shock, standard_atmosphere

Compatibility re-exports in older modules are listed by their defining module
rather than documented a second time.

.. automodule:: aerophysics

.. toctree::
   :maxdepth: 1

   thermophysical
   compressible_flow
   viscous_flow
   flight_units
   errors
   ../references
