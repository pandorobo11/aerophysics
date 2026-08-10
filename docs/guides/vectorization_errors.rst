Vectorization, warnings, and errors
===================================

Vectorized calculations
-----------------------

Public calculation functions accept NumPy-compatible array-like inputs and
broadcast compatible shapes:

>>> import numpy as np
>>> from aerophysics import standard_atmosphere
>>> state = standard_atmosphere(np.array([[0.0], [10_000.0]]))
>>> state.pressure.shape
(2, 1)
>>> np.round(state.pressure[:, 0], 1)
array([101325. ,  26499.9])

Scalar calls return Python floats, while array-like calls return float64
arrays:

>>> isinstance(standard_atmosphere(0.0).pressure, float)
True
>>> state.pressure.dtype == np.dtype(np.float64)
True

Range errors
------------

Values outside a model's supported domain raise
:class:`~aerophysics.exceptions.ModelRangeError` when no supported
extrapolation exists:

>>> from aerophysics.exceptions import ModelRangeError
>>> try:
...     standard_atmosphere(100_000.0)
... except ModelRangeError as error:
...     "altitude" in str(error).lower()
True

Applicability warnings
----------------------

Some correlations remain numerically evaluable outside a fitted range and
emit :class:`~aerophysics.exceptions.ApplicabilityWarning`. Treat the warning
as a request for independent evidence, not as confirmation that the result is
physically valid:

>>> import warnings
>>> from aerophysics.exceptions import ApplicabilityWarning
>>> from aerophysics.transport import AIR_KEYES_VISCOSITY
>>> with warnings.catch_warnings(record=True) as caught:
...     warnings.simplefilter("always")
...     _ = AIR_KEYES_VISCOSITY.dynamic_viscosity(2000.0)
>>> any(issubclass(item.category, ApplicabilityWarning) for item in caught)
True

An oblique-shock request with no attached solution instead raises
:class:`~aerophysics.exceptions.NoAttachedShockError`; it never silently
switches to a normal shock. See :doc:`../models/shock_waves` and
:doc:`../api/errors` for the complete error types.
