Warnings and exceptions
=======================

.. automodule:: aerophysics.exceptions
   :members:

Use :class:`~aerophysics.exceptions.ModelRangeError` when a model cannot
evaluate the requested state. Treat
:class:`~aerophysics.exceptions.ApplicabilityWarning` as an explicit notice
that a correlation was evaluated beyond its documented evidence. An oblique
shock with no attached solution raises
:class:`~aerophysics.exceptions.NoAttachedShockError` rather than silently
changing the requested physics.

Thermally perfect conical numerical integration, maximization, or root
resolution failures raise :class:`~aerophysics.exceptions.ShockConvergenceError`,
a subclass of ``RuntimeError``. This does not imply temperature-range failure
or physical detachment. The conical GUI displays this failure for a single
calculation; a sweep retains that point with ``status="error"`` and continues
to the remaining points. Unrelated runtime errors are not suppressed.
