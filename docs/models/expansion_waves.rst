Prandtl--Meyer expansion waves
==============================

This model assumes steady, inviscid, adiabatic expansion of a calorically
perfect gas. State ratios are downstream over upstream, and angles are in
radians. Convert explicitly with
:func:`aerophysics.units.degrees_to_radians` and
:func:`aerophysics.units.radians_to_degrees`.

For :math:`M\ge1`, the Prandtl--Meyer function is

.. math::

   \nu(M)
   =\sqrt{\frac{\gamma+1}{\gamma-1}}
    \tan^{-1}\left[
      \sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}
    \right]
    -\tan^{-1}\left(\sqrt{M^2-1}\right).

Its finite limiting value is

.. math::

   \nu_\max=\frac{\pi}{2}
   \left(\sqrt{\frac{\gamma+1}{\gamma-1}}-1\right).

For a centered expansion through turn angle :math:`\delta`,

.. math::

   \nu(M_2)=\nu(M_1)+\delta.

:func:`aerophysics.expansion.mach_from_prandtl_meyer` solves this equation
numerically. Total temperature and total pressure remain constant, while

.. math::

   \frac{T_2}{T_1}=\frac{F(M_1)}{F(M_2)},
   \qquad
   \frac{p_2}{p_1}
   =\left(\frac{T_2}{T_1}\right)^{\gamma/(\gamma-1)},
   \qquad
   \frac{\rho_2}{\rho_1}
   =\left(\frac{T_2}{T_1}\right)^{1/(\gamma-1)}.

.. list-table:: Expansion symbols
   :header-rows: 1
   :widths: 16 29 39 16

   * - Symbol
     - API name
     - Meaning
     - Unit
   * - :math:`\nu`
     - ``prandtl_meyer_angle``
     - Prandtl--Meyer angle
     - rad
   * - :math:`\delta`
     - ``turn_angle``
     - Flow turning angle
     - rad
   * - :math:`M_1,M_2`
     - ``upstream_mach``, ``downstream_mach``
     - Upstream and downstream Mach numbers
     - dimensionless

The input turn must be nonnegative and keep the downstream angle below
:math:`\nu_\max`. The complete state change is returned by
:func:`aerophysics.expansion.prandtl_meyer_expansion`. These governing
relations and reference values follow
:ref:`NACA Report 1135 <ref-naca-report-1135>`.

>>> from aerophysics import prandtl_meyer_expansion
>>> from aerophysics.units import degrees_to_radians
>>> expansion = prandtl_meyer_expansion(2.0, degrees_to_radians(10.0))
>>> round(expansion.downstream_mach, 6)
2.384887
