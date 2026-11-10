Arrays, warnings, and errors
============================

Many cases at once
------------------

Most functions accept arrays as well as single numbers, so you can calculate a
whole sweep in one call. The result has the same shape as the input:

>>> import numpy as np
>>> from aerophysics import standard_atmosphere
>>> state = standard_atmosphere(np.array([[0.0], [10_000.0]]))
>>> state.pressure.shape
(2, 1)
>>> np.round(state.pressure[:, 0], 1)
array([101325. ,  26499.9])

When several inputs are arrays, they are combined by the usual NumPy
broadcasting rules. A single number in gives a Python ``float`` out; an array
in gives a NumPy array out:

>>> isinstance(standard_atmosphere(0.0).pressure, float)
True
>>> state.pressure.dtype == np.dtype(np.float64)
True

Result arrays are read-only, so that the values in one result always belong
together. Make a copy if you need to change one:

>>> state.pressure.flags.writeable
False
>>> writable_pressure = state.pressure.copy()
>>> writable_pressure.flags.writeable
True

What can go wrong
-----------------

The package reports problems in three ways:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Report
     - Meaning
   * - ``ValueError``
     - An input is not usable, for example a negative temperature or a value
       that is not a number. The whole call fails, even if only one element of
       an array is bad.
   * - :class:`~aerophysics.exceptions.ModelRangeError`
     - The inputs make sense, but the model cannot calculate that case, for
       example an altitude above the top of the standard atmosphere.
   * - :class:`~aerophysics.exceptions.ApplicabilityWarning`
     - A result was calculated, but outside the range the formula was fitted
       or checked for. Treat the result with caution.

A few calculations have their own, more specific errors:

.. list-table::
   :header-rows: 1
   :widths: 44 56

   * - Error
     - Meaning
   * - :class:`~aerophysics.exceptions.NoAttachedShockError`
     - The wedge or cone angle is too large for an attached shock.
   * - :class:`~aerophysics.exceptions.BoundaryLayerSeparationError`
     - A marched boundary layer separated. The error holds the results
       upstream of separation; see :doc:`marching_boundary_layers`.
   * - :class:`~aerophysics.exceptions.ShockConvergenceError`,
       :class:`~aerophysics.exceptions.ExpansionConvergenceError`,
       :class:`~aerophysics.exceptions.BoundaryLayerConvergenceError`
     - The numerical solution did not converge. This is rare and does not by
       itself mean that the inputs are wrong.

Range errors
------------

A value outside the range a model supports raises
:class:`~aerophysics.exceptions.ModelRangeError`:

>>> from aerophysics.exceptions import ModelRangeError
>>> try:
...     standard_atmosphere(100_000.0)
... except ModelRangeError as error:
...     "altitude" in str(error).lower()
True

Applicability warnings
----------------------

Some formulas can still be evaluated outside the range they were fitted to.
They then return a value and emit
:class:`~aerophysics.exceptions.ApplicabilityWarning`. The warning means that
you should check the result against another source before relying on it. This
example collects the warnings from a call:

>>> import warnings
>>> from aerophysics.exceptions import ApplicabilityWarning
>>> from aerophysics.transport import AIR_KEYES_VISCOSITY
>>> with warnings.catch_warnings(record=True) as caught:
...     warnings.simplefilter("always")
...     _ = AIR_KEYES_VISCOSITY.dynamic_viscosity(2000.0)
>>> any(issubclass(item.category, ApplicabilityWarning) for item in caught)
True

Python shows each warning once by default. To stop on any applicability
warning, for example in an automated script, turn it into an error::

    import warnings
    from aerophysics.exceptions import ApplicabilityWarning

    warnings.simplefilter("error", ApplicabilityWarning)

Shocks that cannot stay attached
--------------------------------

An oblique or conical shock beyond its maximum angle raises
:class:`~aerophysics.exceptions.NoAttachedShockError`. The package never
replaces it with a normal shock or another model without telling you:

>>> from math import radians
>>> from aerophysics import oblique_shock
>>> from aerophysics.exceptions import NoAttachedShockError
>>> try:
...     oblique_shock(2.0, radians(30.0))
... except NoAttachedShockError:
...     print("no attached shock")
no attached shock

All error types are listed in :doc:`../api/errors`.
