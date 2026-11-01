Warnings and exceptions
=======================

.. automodule:: aerophysics.exceptions
   :members:

The package separates three kinds of problems:

* :class:`~aerophysics.exceptions.ModelRangeError`: the model cannot evaluate
  the requested state, for example because a temperature lies outside the
  gas data.
* :class:`~aerophysics.exceptions.ApplicabilityWarning`: a correlation was
  evaluated, but outside the range supported by its source. Check the result
  against independent evidence.
* :class:`~aerophysics.exceptions.NoAttachedShockError`: the requested wedge
  or cone angle is larger than the attached-shock limit. The request is not
  silently changed into a normal shock or another model.

Numerical failures of the iterative solvers raise
:class:`~aerophysics.exceptions.ShockConvergenceError` or
:class:`~aerophysics.exceptions.ExpansionConvergenceError`. Both are
subclasses of ``RuntimeError``. They mean that the solver could not reach a
reliable answer; they do not mean that the state is out of range or that the
shock has detached. See :doc:`../guides/vectorization_errors` for examples.
