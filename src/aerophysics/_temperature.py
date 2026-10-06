"""Roundoff repair only for temperatures restored from validated solver ratios."""

import numpy as np

from aerophysics.gas import PerfectGas
from aerophysics.real_gas import BeattieBridgemanGas, HarmonicOscillatorGas
from aerophysics.thermochemistry import ThermallyPerfectGas


def restore_static_temperature(
    temperature: float,
    gas: PerfectGas | ThermallyPerfectGas | HarmonicOscillatorGas | BeattieBridgemanGas,
    *,
    allow_extrapolation: bool = False,
) -> float:
    """Undo one-ULP range overshoot when restoring a solved T from a ratio."""
    if not allow_extrapolation and isinstance(
        gas, (ThermallyPerfectGas, HarmonicOscillatorGas)
    ):
        minimum, maximum = gas.temperature_range
        if maximum < temperature <= np.nextafter(maximum, np.inf):
            return maximum
        if np.nextafter(minimum, -np.inf) <= temperature < minimum:
            return minimum
    return temperature
