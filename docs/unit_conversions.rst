Unit conversions
================

The functions in :mod:`aerophysics.units` make customary-to-SI conversions
explicit at the package boundary. Aerodynamic calculation APIs otherwise use
SI units. Forward and inverse functions share the same exact constants.

Conversion equations
--------------------

.. list-table:: Implemented forward conversions
   :header-rows: 1
   :widths: 34 44 22

   * - Conversion
     - Equation
     - API
   * - International feet to metres
     - :math:`x_\mathrm{m}=0.3048x_\mathrm{ft}`
     - ``feet_to_meters``
   * - Knots to metres per second
     - :math:`V_\mathrm{m/s}=(1852/3600)V_\mathrm{kt}`
     - ``knots_to_meters_per_second``
   * - Degrees Fahrenheit to kelvin
     - :math:`T_\mathrm{K}=(T_\mathrm{^\circ F}-32)(5/9)+273.15`
     - ``fahrenheit_to_kelvin``
   * - Pounds-force per square inch to pascals
     - :math:`p_\mathrm{Pa}=6894.757293168p_\mathrm{psi}`
     - ``psi_to_pascals``
   * - Pounds-force per square foot to pascals
     - :math:`p_\mathrm{Pa}=47.8802589803p_\mathrm{psf}`
     - ``psf_to_pascals``
   * - Pounds mass to kilograms
     - :math:`m_\mathrm{kg}=0.45359237m_\mathrm{lbm}`
     - ``pounds_mass_to_kilograms``
   * - Slugs to kilograms
     - :math:`m_\mathrm{kg}=14.5939029372m_\mathrm{slug}`
     - ``slugs_to_kilograms``
   * - Degrees to radians
     - :math:`\theta_\mathrm{rad}=(\pi/180)\theta_\mathrm{deg}`
     - ``degrees_to_radians``

Inverse conversions divide by the corresponding scale factor. Temperature
uses the affine inverse

.. math::

   T_\mathrm{^\circ F}=(T_\mathrm{K}-273.15)(9/5)+32.

Values below absolute zero are rejected by the temperature converters.

Examples
--------

>>> from aerophysics.units import feet_to_meters
>>> from aerophysics.units import knots_to_meters_per_second
>>> from aerophysics.units import degrees_to_radians
>>> feet_to_meters(10_000.0)
3048.0
>>> round(knots_to_meters_per_second(100.0), 6)
51.444444
>>> round(degrees_to_radians(90.0), 8)
1.57079633

The same functions accept NumPy-compatible array-like values. See the
:mod:`aerophysics.units` API reference for all inverse function names and
return-type behavior.
