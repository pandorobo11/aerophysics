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
