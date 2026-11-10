API reference
=============

The API reference lists every public function and class with its inputs,
outputs, and units, grouped by subject. For worked examples, see the
:doc:`../guides/index`.

The most frequently used classes and functions can be imported directly from
``aerophysics``:

.. code-block:: python

   from aerophysics import FlightCondition, oblique_shock, standard_atmosphere

Everything else is imported from its own module, as shown on each page.

.. automodule:: aerophysics

.. toctree::
   :maxdepth: 1

   thermophysical
   compressible_flow
   viscous_flow
   flight_units
   errors
