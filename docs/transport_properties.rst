Transport properties
====================

The :mod:`aerophysics.transport` module provides interchangeable,
temperature-dependent dynamic-viscosity models. All models accept temperature
in kelvin and return dynamic viscosity in Pa s. Scalar input returns a Python
``float`` and array-like input returns a float64 NumPy array.

Dynamic-viscosity interface
---------------------------

:class:`~aerophysics.transport.DynamicViscosityModel` is the structural
interface consumed by the boundary-layer APIs. Any object implementing
``dynamic_viscosity(temperature)`` with the documented scalar and array
behaviour can be supplied as ``viscosity_model``.

Sutherland model
----------------

:class:`~aerophysics.transport.SutherlandModel` evaluates

.. math::

   \mu(T)=\mu_\mathrm{ref}
   \left(\frac{T}{T_\mathrm{ref}}\right)^{3/2}
   \frac{T_\mathrm{ref}+S}{T+S}.

``AIR_VISCOSITY`` is the existing dry-air default, with
:math:`\mu_\mathrm{ref}=1.7894\times10^{-5}\ \mathrm{Pa\,s}`,
:math:`T_\mathrm{ref}=288.15\ \mathrm{K}`, and :math:`S=110.4\ \mathrm{K}`.
It remains the model used by the standard atmosphere and flight-condition
calculations.

Keyes model
-----------

:class:`~aerophysics.transport.KeyesModel` uses

.. math::

   \mu(T)=\frac{a_0T^{3/2}}
   {T+a_1 10^{-a_2/T}}.

``AIR_KEYES_VISCOSITY`` uses :math:`a_0=1.488\times10^{-6}`,
:math:`a_1=122.1\ \mathrm{K}`, and :math:`a_2=5\ \mathrm{K}`. Its nominal
temperature range is 79--1845 K.

>>> from aerophysics.transport import AIR_KEYES_VISCOSITY
>>> f"{AIR_KEYES_VISCOSITY.dynamic_viscosity(300.0):.7e}"
'1.8519327e-05'

Blottner species and mixture models
-----------------------------------

:class:`~aerophysics.transport.BlottnerModel` represents one gas species:

.. math::

   \mu_s(T)=0.1\exp\left[(A_s\ln T+B_s)\ln T+C_s\right].

The nominal fitted range is 1000--30000 K. The built-in frozen dry-air model
uses the following coefficients.

.. list-table:: Built-in Blottner coefficients
   :header-rows: 1

   * - Species
     - :math:`A_s`
     - :math:`B_s`
     - :math:`C_s`
   * - N₂
     - 0.0268142
     - 0.3177838
     - -11.3155513
   * - O₂
     - 0.0449290
     - -0.0826158
     - -9.2019475
   * - Ar
     - -0.02201
     - 1.010
     - -13.42
   * - CO₂
     - -0.041372
     - 1.3293
     - -15.016

:class:`~aerophysics.transport.WilkeMixtureViscosityModel` combines species
viscosities using mole fractions :math:`x_i`, molar masses :math:`M_i`, and
Wilke's rule:

.. math::

   \mu=\sum_i\frac{x_i\mu_i}{\sum_jx_j\phi_{ij}},
   \qquad
   \phi_{ij}=\frac{\left[1+\left(\frac{\mu_i}{\mu_j}\right)^{1/2}
   \left(\frac{M_j}{M_i}\right)^{1/4}\right]^2}
   {\sqrt{8\left(1+M_i/M_j\right)}}.

``AIR_BLOTTNER_VISCOSITY`` uses the same normalized N₂/O₂/Ar/CO₂ mole
fractions and molar masses as ``AIR_NASA7`` and ``AIR_NASA9``.

>>> from aerophysics.transport import AIR_BLOTTNER_VISCOSITY
>>> f"{AIR_BLOTTNER_VISCOSITY.dynamic_viscosity(1000.0):.7e}"
'4.1375747e-05'

Temperatures outside a Keyes or Blottner nominal range remain computable but
emit :class:`~aerophysics.exceptions.ApplicabilityWarning`. Non-positive or
non-finite temperatures raise ``ValueError``. The Blottner dry-air preset has
fixed composition: it does not model dissociation, reactions, ionization, or
thermodynamic nonequilibrium.

Thermal conductivity
--------------------

:class:`~aerophysics.transport.USSAConductivityModel` evaluates the U.S.
Standard Atmosphere correlation

.. math::

   k(T)=\frac{c_kT^{3/2}}{T+A_k10^{-B_k/T}},

with :math:`c_k=2.64638\times10^{-3}`, :math:`A_k=245.4\ \mathrm{K}`, and
:math:`B_k=12\ \mathrm{K}` for ``AIR_CONDUCTIVITY``. The result is in W/(m K).

The model equations, coefficients, and source reports are listed in
:doc:`references`.
