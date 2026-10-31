"""Analytical and regression checks, not independent experimental validation."""

import warnings
from dataclasses import FrozenInstanceError, fields
from typing import Any

import numpy as np
import pytest

from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.gas import AIR, PerfectGas
from aerophysics.heating import (
    DEFAULT_ATTACHMENT_LINE_RECOVERY_FACTOR,
    POLL_TURBULENT_COEFFICIENT,
    POLL_TURBULENT_EXPONENT,
    AttachmentLineHeatingModel,
    AttachmentLineHeatingResult,
    _turbulent_attachment_line,
)
from aerophysics.heating import (
    swept_cylinder_turbulent_heat_flux as cylinder,
)
from aerophysics.transport import AIR_VISCOSITY, SutherlandModel

ARGS = (10.0, 226.5, 1197.0, 0.1, np.pi / 6)


def evaluate(*args: Any, **kwargs: Any) -> AttachmentLineHeatingResult:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        return cylinder(*(args or ARGS), **kwargs)


@pytest.mark.parametrize("angle", np.deg2rad([15, 30, 60]))
def test_cr_181893_equation_512(angle: float) -> None:
    result = evaluate(*ARGS[:4], angle, wall_temperature=300)
    diameter = 2 * ARGS[3]
    velocity = result.freestream_velocity
    rho_inf = result.freestream_density
    reference = result.reference_temperature
    mu_inf = AIR_VISCOSITY.dynamic_viscosity(ARGS[1])
    mu_ref = AIR_VISCOSITY.dynamic_viscosity(reference)
    re_diameter = rho_inf * velocity * diameter / mu_inf
    u1_zero_sweep = result.velocity_gradient * diameter / (velocity * np.cos(angle))
    # Independent freestream normalization, including diameter and rho* basis.
    st_inf = (
        0.0689
        / (2 * 0.72 ** (2 / 3))
        * (result.edge_density * result.edge_temperature / (rho_inf * reference))
        ** 0.79
        * (mu_ref * u1_zero_sweep / (mu_inf * re_diameter)) ** 0.21
        * np.sin(angle) ** 0.58
        * np.cos(angle) ** 0.21
    )
    flux = st_inf * rho_inf * velocity * AIR.cp * (result.recovery_temperature - 300)
    assert result.heat_flux == pytest.approx(flux, rel=1e-12)


@pytest.mark.parametrize("recovery_factor", [0.89, 1.2])
def test_core_hand_calculation(recovery_factor: float) -> None:
    te, pe, we, gradient, wall = (
        np.asarray(x) for x in (900.0, 50000.0, 1200.0, 8000.0, 350.0)
    )
    # Custom gas and transport also verify that the supplied models are used.
    gas = PerfectGas(heat_capacity_ratio=1.3, specific_gas_constant=300)
    viscosity = SutherlandModel(1.8e-5, 300, 110)
    actual = _turbulent_attachment_line(
        te,
        pe,
        we,
        gradient,
        wall,
        prandtl_number=0.8,
        recovery_factor=recovery_factor,
        gas=gas,
        viscosity_model=viscosity,
    )
    taw = te + recovery_factor * we**2 / (2 * gas.cp)
    reference = te + 0.1 * (wall - te) + 0.6 * (taw - te)
    rho = pe / (gas.specific_gas_constant * reference)
    mu = 1.8e-5 * (reference / 300) ** 1.5 * (300 + 110) / (reference + 110)
    reynolds = we / np.sqrt(mu / rho * gradient)
    stanton = 0.0689 / (2 * 0.8 ** (2 / 3)) * reynolds**-0.42
    flux = stanton * rho * we * gas.cp * (taw - wall)
    me = we / np.sqrt(1.3 * 300 * te)
    assert actual == pytest.approx(
        (flux, taw, reference, reynolds, stanton, me), rel=1e-12
    )


def test_edge_state_and_energy() -> None:
    result = evaluate(wall_temperature=300)
    m, t, p, radius, sweep = ARGS
    gamma = AIR.heat_capacity_ratio
    normal_mach = m * np.cos(sweep)
    # Independent normal-shock p02/p01 relation; using full M fails this test.
    pressure_ratio = (
        (gamma + 1) * normal_mach**2 / ((gamma - 1) * normal_mach**2 + 2)
    ) ** (gamma / (gamma - 1)) * (
        (gamma + 1) / (2 * gamma * normal_mach**2 - (gamma - 1))
    ) ** (1 / (gamma - 1))
    factor = 1 + (gamma - 1) / 2 * normal_mach**2
    assert result.edge_temperature == pytest.approx(t * factor, rel=1e-12)
    assert result.edge_pressure == pytest.approx(
        p * factor ** (gamma / (gamma - 1)) * pressure_ratio, rel=1e-12
    )
    assert result.edge_density == pytest.approx(
        result.edge_pressure / (AIR.specific_gas_constant * result.edge_temperature),
        rel=1e-12,
    )
    assert result.edge_velocity == pytest.approx(
        result.freestream_velocity * np.sin(sweep), rel=1e-12
    )
    assert result.velocity_gradient == pytest.approx(
        np.sqrt(2 * (result.edge_pressure - p) / result.edge_density) / radius,
        rel=1e-12,
    )
    assert (
        AIR.cp * result.edge_temperature + result.edge_velocity** 2 / 2
        == pytest.approx(AIR.cp * t + result.freestream_velocity**2 / 2, rel=1e-12)
    )
    assert result.edge_mach == pytest.approx(
        result.edge_velocity
        / np.sqrt(gamma * AIR.specific_gas_constant * result.edge_temperature),
        rel=1e-12,
    )


def test_sign_scaling_and_normalization() -> None:
    base = evaluate(wall_temperature=300)
    assert evaluate(
        wall_temperature=base.recovery_temperature
    ).heat_flux == pytest.approx(0, abs=1e-8)
    assert evaluate(wall_temperature=base.recovery_temperature + 100).heat_flux < 0
    large = evaluate(*ARGS[:3], 0.4, ARGS[4], wall_temperature=300)
    assert large.velocity_gradient == pytest.approx(
        base.velocity_gradient / 4, rel=1e-12
    )
    assert large.heat_flux == pytest.approx(base.heat_flux * 4**-0.21, rel=1e-12)
    rho = base.edge_pressure / (AIR.specific_gas_constant * base.reference_temperature)
    nu = AIR_VISCOSITY.dynamic_viscosity(base.reference_temperature) / rho
    assert base.reynolds_number == pytest.approx(
        base.edge_velocity / np.sqrt(nu * base.velocity_gradient), rel=1e-12
    )
    assert base.heat_flux == pytest.approx(
        base.stanton_number
        * rho
        * base.edge_velocity
        * AIR.cp
        * (base.recovery_temperature - 300),
        rel=1e-12,
    )


def test_regression_values() -> None:
    result = evaluate(wall_temperature=300)
    # Regression values only; these are not experimental validation.
    expected = {
        "heat_flux": 2287189.537895632,
        "reynolds_number": 431.35648833491416,
        "stanton_number": 0.00335482317626368,
        "edge_temperature": 3624.0,
        "edge_pressure": 116142.94709632447,
        "velocity_gradient": 14349.626672688486,
    }
    for name, value in expected.items():
        assert getattr(result, name) == pytest.approx(value, rel=1e-10)
    assert result.model is AttachmentLineHeatingModel.POLL
    assert result.sweep_angle == ARGS[4]
    assert DEFAULT_ATTACHMENT_LINE_RECOVERY_FACTOR == 0.89
    assert POLL_TURBULENT_COEFFICIENT == 0.0689
    assert POLL_TURBULENT_EXPONENT == 0.42


def test_arrays_and_immutability() -> None:
    mach = np.array([[6.0], [7.0]])
    result = evaluate(mach, 226.5, 1197, [0.1, 0.4], np.pi / 6, wall_temperature=300)
    for i, m in enumerate((6, 7)):
        for j, radius in enumerate((0.1, 0.4)):
            scalar = evaluate(m, 226.5, 1197, radius, np.pi / 6, wall_temperature=300)
            for field in fields(result):
                array_value = getattr(result, field.name)
                scalar_value = getattr(scalar, field.name)
                if field.name == "model":
                    assert (
                        array_value is scalar_value is AttachmentLineHeatingModel.POLL
                    )
                else:
                    assert isinstance(scalar_value, float)
                    assert isinstance(array_value, np.ndarray)
                    assert array_value.shape == (2, 2)
                    assert array_value.dtype == np.float64
                    assert not array_value.flags.writeable
                    assert array_value[i, j] == pytest.approx(scalar_value, rel=1e-12)
    mach[:] = 99
    assert result.freestream_velocity[0, 0] == pytest.approx(1810.2162990262989)  # type: ignore[index]
    with pytest.raises(FrozenInstanceError):
        result.heat_flux = 0  # type: ignore[misc]


@pytest.mark.parametrize("sweep", [0, np.pi / 2, -0.1, np.pi / 2 + 0.1, np.nan, np.inf])
def test_invalid_sweep(sweep: float) -> None:
    with pytest.raises(ValueError, match="sweep_angle"):
        cylinder(*ARGS[:4], sweep, wall_temperature=300)


def test_invalid_normal_mach() -> None:
    for mach, angle in [(2, 1.1), (2, np.pi / 3 + 1e-14), (1 / np.cos(0.5), 0.5)]:
        with pytest.raises(ValueError, match="normal Mach number"):
            cylinder(mach, *ARGS[1:4], angle, wall_temperature=300)


@pytest.mark.parametrize("index", [0, 1, 2, 3, 5])
@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
def test_invalid_common_inputs(index: int, value: float) -> None:
    args: Any = [*ARGS, 300]
    args[index] = value
    with pytest.raises(ValueError):
        cylinder(*args[:5], wall_temperature=args[5])


@pytest.mark.parametrize("name", ["prandtl_number", "recovery_factor"])
@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
def test_invalid_parameters(name: str, value: float) -> None:
    kwargs: Any = {name: value}
    with pytest.raises(ValueError, match=name):
        cylinder(*ARGS, wall_temperature=300, **kwargs)


def test_other_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="mach must be greater than one"):
        cylinder(1, *ARGS[1:], wall_temperature=300)
    with pytest.raises(ValueError, match="model must be an AttachmentLineHeatingModel"):
        cylinder(*ARGS, wall_temperature=300, model="poll")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="heating inputs must be broadcastable"):
        cylinder([6, 7], [200, 250, 300], *ARGS[2:], wall_temperature=300)


@pytest.mark.parametrize("array", [False, True])
@pytest.mark.parametrize(
    "mach,sweep,wall,expected",
    [
        (10.0, np.pi / 6, 1000.0, ["2000 K"]),
        (6.0, 0.05, 300.0, ["edge Mach number"]),
        (20.0, 1.5, 30000.0, ["2000 K", "edge Mach number", "wall-to-total"]),
        (6.0, np.pi / 6, 100.0, ["wall-to-total"]),
        (6.0, np.pi / 6, 1900.0, ["wall-to-total"]),
        (6.0, np.pi / 6, 300.0, []),
        # A hot wall within the database has no extra warning.
        (6.0, np.pi / 3, 1750.0, []),
    ],
)
def test_applicability_warnings(
    mach: float, sweep: float, wall: float, expected: list[str], array: bool
) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cylinder(
            [mach, mach] if array else mach,
            226.5,
            1197,
            0.1,
            sweep,
            wall_temperature=wall,
        )
    assert len(caught) == len(expected)
    for warning, message in zip(caught, expected, strict=True):
        assert warning.category is ApplicabilityWarning
        assert message in str(warning.message)


@pytest.mark.parametrize("edge_mach", [0.35, 6.89])
@pytest.mark.parametrize("wall_ratio", [0.14, 0.96])
def test_closed_database_warning_boundaries(
    monkeypatch: pytest.MonkeyPatch, edge_mach: float, wall_ratio: float
) -> None:
    # Isolate the warning boundary from edge-state floating-point inversion.
    def core(*args: Any, **kwargs: Any) -> Any:
        values = _turbulent_attachment_line(*args, **kwargs)
        return (*values[:-1], np.asarray(edge_mach))

    monkeypatch.setattr("aerophysics.heating._turbulent_attachment_line", core)
    factor = 1 + (AIR.heat_capacity_ratio - 1) / 2 * 6**2
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cylinder(
            [6, 6],
            1000 / factor,
            1197,
            0.1,
            np.pi / 6,
            wall_temperature=wall_ratio * 1000,
        )
    assert not caught


def test_total_temperature_warning_boundary() -> None:
    factor = 1 + (AIR.heat_capacity_ratio - 1) / 2 * 6**2
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cylinder(6, 2000 / factor, 1197, 0.1, np.pi / 6, wall_temperature=500)
    assert not caught
