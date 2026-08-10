Thermally perfect gases
========================

``aerophysics`` provides NASA seven- and nine-coefficient polynomial
evaluators and a frozen-composition ideal-gas mixture model. These models make
heat capacity a function of temperature while retaining the ideal-gas equation
of state. They do not add chemical equilibrium or thermodynamic
nonequilibrium.

The seven-coefficient form and data are documented by
:ref:`McBride, Gordon, and Reno (1993) <ref-mcbride-gordon-reno-1993>`; the
nine-coefficient form is documented by :ref:`McBride, Zehe, and Gordon (2002)
<ref-mcbride-zehe-gordon-2002>`.

NASA polynomial forms
---------------------

For the NASA seven-coefficient form, each temperature region stores
:math:`a_1,\ldots,a_7` and evaluates

.. math::

   \frac{\bar c_p^\circ}{R_u}
   &= a_1+a_2T+a_3T^2+a_4T^3+a_5T^4,\\
   \frac{\bar h^\circ}{R_uT}
   &= a_1+\frac{a_2T}{2}+\frac{a_3T^2}{3}
      +\frac{a_4T^3}{4}+\frac{a_5T^4}{5}+\frac{a_6}{T},\\
   \frac{\bar s^\circ}{R_u}
   &= a_1\ln T+a_2T+\frac{a_3T^2}{2}
      +\frac{a_4T^3}{3}+\frac{a_5T^4}{4}+a_7.

The NASA nine-coefficient form adds inverse-temperature terms:

.. math::

   \frac{\bar c_p^\circ}{R_u}
   &= a_1T^{-2}+a_2T^{-1}+a_3+a_4T+a_5T^2+a_6T^3+a_7T^4,\\
   \frac{\bar h^\circ}{R_uT}
   &= -a_1T^{-2}+a_2\frac{\ln T}{T}+a_3+\frac{a_4T}{2}
      +\frac{a_5T^2}{3}+\frac{a_6T^3}{4}
      +\frac{a_7T^4}{5}+\frac{a_8}{T},\\
   \frac{\bar s^\circ}{R_u}
   &= -\frac{a_1T^{-2}}{2}-a_2T^{-1}+a_3\ln T+a_4T
      +\frac{a_5T^2}{2}+\frac{a_6T^3}{3}
      +\frac{a_7T^4}{4}+a_9.

``NASA7Polynomial`` and ``NASA9Polynomial`` return these three dimensionless
quantities through ``cp_over_r``, ``h_over_rt``, and ``s_over_r``. A shared
temperature boundary belongs to the lower-temperature region. Inputs outside
the fitted range raise :class:`~aerophysics.exceptions.ModelRangeError`.
Passing ``allow_extrapolation=True`` uses the nearest end region and emits
:class:`~aerophysics.exceptions.ApplicabilityWarning`.

Species and mixture properties
------------------------------

An :class:`~aerophysics.thermochemistry.IdealGasSpecies` combines a polynomial
with its molar mass and standard-state pressure. Species methods return molar
SI properties. For fixed mole fractions :math:`x_i`, the mixture relations are

.. math::

   M &= \sum_i x_i M_i,
   &R &= \frac{R_u}{M},\\
   c_p &= \frac{\sum_i x_i\bar c_{p,i}^\circ}{M},
   &c_v &= c_p-R,\\
   \gamma &= \frac{c_p}{c_v},
   &a &= \sqrt{\gamma RT}.

The ideal-mixture entropy at pressure :math:`p` is

.. math::

   s(T,p)=\frac{1}{M}\sum_i x_i
   \left[
      \bar s_i^\circ(T)
      -R_u\ln\left(\frac{x_i p}{p_i^\circ}\right)
   \right].

NASA standard enthalpy includes formation enthalpy. It is therefore distinct
from the sensible enthalpy used for a temperature change. The
``standard_enthalpy`` and ``standard_internal_energy`` methods preserve the
NASA reference, while ``sensible_enthalpy`` and
``sensible_internal_energy`` subtract the value at an explicit reference
temperature.

Built-in frozen dry air
-----------------------

``AIR_NASA7`` and ``AIR_NASA9`` use normalized U.S. Standard Atmosphere dry-air
mole fractions for N2, O2, Ar, and CO2. Both coefficient sets have a common
fitted range of 200--6000 K:

``AIR_NASA7`` uses the :ref:`Cantera NASA gas data
<ref-cantera-nasa-gas-data>` and ``AIR_NASA9`` uses the
:ref:`NASA CEA data <ref-nasa-cea-data>` pinned by the package.

>>> from aerophysics import AIR_NASA9
>>> round(AIR_NASA9.cp(300.0), 3)
1004.829
>>> round(AIR_NASA9.cv(1000.0), 3)
853.978
>>> round(AIR_NASA9.heat_capacity_ratio(2000.0), 6)
1.297997
>>> round(AIR_NASA9.sensible_enthalpy(1000.0), 1)
747891.1

The two models are useful for comparing polynomial parameterizations. The
isentropic API accepts them together with total temperature and integrates
their enthalpy and entropy variation. Shock and Prandtl--Meyer APIs remain
calorically perfect. Substituting a local :math:`\gamma(T)` into a
constant-:math:`\gamma` relation does not make that relation thermally
perfect.

At high temperature, real air dissociates and eventually ionizes. The built-in
models deliberately keep N2, O2, Ar, and CO2 mole fractions fixed through
6000 K, so their high-temperature results describe a frozen-composition
ideal-gas model rather than equilibrium air.
