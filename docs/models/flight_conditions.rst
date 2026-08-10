Flight conditions
=================

:class:`aerophysics.flight.FlightCondition` combines a standard-atmosphere
state with either a Mach number or a velocity. It is the convenient entry
point when one consistent ambient, dynamic, Reynolds-number, and total state
is needed. Atmospheric properties follow :ref:`U.S. Standard Atmosphere 1976
<ref-us-standard-atmosphere-1976>`, and the calorically perfect total-state
relations follow :ref:`NACA Report 1135 <ref-naca-report-1135>`.

Relations
---------

At the selected geometric altitude,

.. math::

   V=Ma,
   \qquad
   M=\frac{V}{a},
   \qquad
   q=\frac{1}{2}\rho V^2.

The Reynolds number per unit length and, when a positive characteristic
length :math:`L` is supplied, the dimensionless Reynolds number are

.. math::

   \frac{Re_L}{L}=\frac{\rho V}{\mu},
   \qquad
   Re_L=\frac{\rho VL}{\mu}.

Define :math:`F(M)=1+(\gamma-1)M^2/2`. The isentropic total conditions are

.. math::

   T_0=TF(M),
   \qquad
   p_0=pF(M)^{\gamma/(\gamma-1)},
   \qquad
   \rho_0=\rho F(M)^{1/(\gamma-1)}.

The atmosphere supplies :math:`T`, :math:`p`, :math:`\rho`, :math:`a`, and
:math:`\mu`; see :doc:`gas_and_atmosphere`. The total-state convention and
factor are described in :doc:`isentropic_flow`.

.. list-table:: Flight-condition symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - SI unit
   * - :math:`M`
     - ``mach``
     - Mach number
     - dimensionless
   * - :math:`V`
     - ``velocity``
     - Flight speed
     - m/s
   * - :math:`q`
     - ``dynamic_pressure``
     - Dynamic pressure
     - Pa
   * - :math:`L`
     - ``characteristic_length``
     - Optional reference length
     - m
   * - :math:`Re_L/L`
     - ``reynolds_number_per_length``
     - Reynolds number per unit length
     - 1/m
   * - :math:`Re_L`
     - ``reynolds_number``
     - Reynolds number based on :math:`L`
     - dimensionless
   * - :math:`T_0,p_0,\rho_0`
     - ``total_temperature``, ``total_pressure``, ``total_density``
     - Isentropic total state
     - K, Pa, kg/m³

Inputs and applicability
------------------------

Use :meth:`~aerophysics.flight.FlightCondition.from_mach` or
:meth:`~aerophysics.flight.FlightCondition.from_velocity`; do not supply two
independent speed definitions. The altitude must lie within the standard
atmosphere range. ``reynolds_number`` is ``None`` unless a positive
``characteristic_length`` is supplied, while ``reynolds_number_per_length``
is always available.

>>> from aerophysics import FlightCondition
>>> condition = FlightCondition.from_mach(
...     0.0, 0.8, characteristic_length=1.5
... )
>>> round(condition.velocity, 3)
272.235
>>> round(condition.dynamic_pressure, 1)
45393.6
>>> round(condition.total_temperature, 4)
325.0332
>>> round(condition.reynolds_number)
27955292
