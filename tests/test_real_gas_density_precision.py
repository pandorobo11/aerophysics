"""Relative-precision checks for the public dilute-gas density inverse."""

import pytest

from aerophysics.real_gas import BeattieBridgemanGas


@pytest.mark.parametrize("pressure", [2.5e-303, 1.0e-200, 1.0e-9, 1.0, 1.0e5, 1.0e7])
def test_dilute_density_inverse_and_eos_closure_use_relative_tolerances(
    pressure: float,
) -> None:
    gas = BeattieBridgemanGas(287.0, 1.4, 0.0, 0.0, 0.0, 0.0, 0.0)
    temperature = 300.0
    density = float(gas.density(temperature, pressure))
    expected_density = pressure / (gas.specific_gas_constant * temperature)

    # pytest.approx defaults to an absolute tolerance of 1e-12, which can
    # otherwise accept a large relative error at these dilute-state scales.
    assert density == pytest.approx(expected_density, rel=1e-14, abs=0.0)
    closure = float(gas.pressure(temperature, density)) / pressure
    assert closure == pytest.approx(1.0, rel=1e-14, abs=0.0)
