"""Physical checks for local state and transport columns in flow calculators."""

from typing import Any

import numpy as np
import pytest

from aerophysics import AIR_HARMONIC_OSCILLATOR, AIR_NASA7, AIR_NASA9
from aerophysics.gas import AIR, PerfectGas
from aerophysics.gui.adapters import (
    CalculationResult,
    Row,
    conical_shock_condition,
    conical_shock_sweep,
    isentropic_condition,
    normal_shock_condition,
    normal_shock_sweep,
    oblique_shock_condition,
    oblique_shock_sweep,
)
from aerophysics.gui.tables import display_rows, rows_to_csv
from aerophysics.gui.units import UnitPreferences
from aerophysics.shocks import ShockBranch, ShockGasModel
from aerophysics.transport import AIR_KEYES_VISCOSITY, AIR_VISCOSITY

GASES: dict[str, ShockGasModel] = {
    "AIR": AIR,
    "NASA7": AIR_NASA7,
    "NASA9": AIR_NASA9,
    "HARMONIC_OSCILLATOR": AIR_HARMONIC_OSCILLATOR,
}


def number(row: Row, key: str) -> float:
    value = row[key]
    assert isinstance(value, float)
    return value


def shock_case(kind: str, **options: Any) -> CalculationResult:
    if kind == "normal":
        return normal_shock_condition(upstream_mach=3.0, **options)
    if kind == "oblique":
        return oblique_shock_condition(
            upstream_mach=3.0,
            deflection_angle=float(np.deg2rad(10.0)),
            branch=ShockBranch.WEAK,
            **options,
        )
    return conical_shock_condition(
        upstream_mach=3.0,
        cone_half_angle=float(np.deg2rad(10.0)),
        **options,
    )


@pytest.mark.parametrize("kind", ["normal", "oblique", "cone"])
@pytest.mark.parametrize("gas_model", list(GASES))
def test_shock_output_states_conserve_energy_entropy_and_mass(
    kind: str, gas_model: str
) -> None:
    row = shock_case(
        kind,
        gas_model=gas_model,
        upstream_temperature=500.0,
        upstream_pressure=25_000.0,
        characteristic_length=0.2,
        with_heat_capacities=True,
    ).rows[0]
    gas = GASES[gas_model]
    end = "surface" if kind == "cone" else "downstream"
    v1, v2 = number(row, "upstream_velocity"), number(row, f"{end}_velocity")
    t1, t2 = number(row, "upstream_temperature"), number(row, f"{end}_temperature")
    p1, p2 = number(row, "upstream_pressure"), number(row, f"{end}_pressure")
    rho1, rho2 = number(row, "upstream_density"), number(row, f"{end}_density")
    assert p2 == pytest.approx(rho2 * gas.specific_gas_constant * t2, rel=2e-12)
    assert number(row, "velocity_ratio") == pytest.approx(v2 / v1, rel=2e-12)
    assert number(row, "dynamic_pressure_ratio") == pytest.approx(
        rho2 * v2**2 / (rho1 * v1**2), rel=2e-12
    )
    for prefix, temperature, density, velocity in [
        ("upstream", t1, rho1, v1),
        (end, t2, rho2, v2),
    ]:
        gamma = number(row, f"{prefix}_heat_capacity_ratio")
        assert gamma == pytest.approx(
            number(row, f"{prefix}_cp") / number(row, f"{prefix}_cv"), rel=2e-12
        )
        sound_speed = float(gas.speed_of_sound(temperature))
        assert number(row, f"{prefix}_speed_of_sound") == pytest.approx(
            sound_speed, rel=2e-12
        )
        mu = float(AIR_VISCOSITY.dynamic_viscosity(temperature))
        assert number(row, f"{prefix}_dynamic_viscosity") == pytest.approx(
            mu, rel=2e-12
        )
        assert number(row, f"{prefix}_reynolds_number_per_length") == pytest.approx(
            density * velocity / mu, rel=2e-12
        )
        assert number(row, f"{prefix}_reynolds_number") == pytest.approx(
            0.2 * density * velocity / mu, rel=2e-12
        )
    if isinstance(gas, PerfectGas):
        delta_h = AIR.cp * (t2 - t1)
        delta_s_over_r = AIR.cp / AIR.specific_gas_constant * np.log(t2 / t1) - np.log(
            p2 / p1
        )
    else:
        delta_h = float(gas.standard_enthalpy(t2)) - float(gas.standard_enthalpy(t1))
        delta_s_over_r = (
            float(gas.entropy(t2, p2)) - float(gas.entropy(t1, p1))
        ) / gas.specific_gas_constant
    assert delta_h + 0.5 * (v2**2 - v1**2) == pytest.approx(0.0, abs=2e-3)
    assert number(row, "entropy_change_over_r") == pytest.approx(
        delta_s_over_r, abs=2e-8
    )
    if kind == "normal":
        assert rho1 * v1 == pytest.approx(rho2 * v2, rel=2e-10)
        assert number(row, "velocity_ratio") == pytest.approx(
            1.0 / number(row, "static_density_ratio"), rel=2e-10
        )
    elif kind == "oblique":
        beta, theta = number(row, "shock_angle"), number(row, "deflection_angle")
        assert rho1 * v1 * np.sin(beta) == pytest.approx(
            rho2 * v2 * np.sin(beta - theta), rel=2e-10
        )
        assert v1 * np.cos(beta) == pytest.approx(v2 * np.cos(beta - theta), rel=2e-10)
        assert number(row, "velocity_ratio") != pytest.approx(
            1.0 / number(row, "static_density_ratio"), rel=1e-3
        )


@pytest.mark.parametrize("kind", ["normal", "oblique", "cone"])
def test_pressure_scales_only_absolute_density_pressure_and_reynolds(kind: str) -> None:
    options: dict[str, Any] = dict(gas_model="NASA9", upstream_temperature=500.0)
    a = shock_case(kind, upstream_pressure=25_000.0, **options).rows[0]
    b = shock_case(kind, upstream_pressure=50_000.0, **options).rows[0]
    no_pressure = shock_case(kind, **options).rows[0]
    end = "surface" if kind == "cone" else "downstream"
    for prefix in ["upstream", end]:
        for field in [
            "pressure",
            "density",
            "dynamic_pressure",
            "reynolds_number_per_length",
        ]:
            assert number(b, f"{prefix}_{field}") == pytest.approx(
                2.0 * number(a, f"{prefix}_{field}"), rel=2e-12
            )
        assert number(b, f"{prefix}_velocity") == number(
            no_pressure, f"{prefix}_velocity"
        )
    assert "upstream_reynolds_number_per_length" not in no_pressure
    assert "upstream_cp" not in no_pressure
    assert b["velocity_ratio"] == no_pressure["velocity_ratio"]


def test_legacy_air_ratios_and_optional_static_temperature() -> None:
    old = normal_shock_condition(upstream_mach=3.0).rows[0]
    assert old["upstream_heat_capacity_ratio"] == 1.4
    assert "upstream_velocity" not in old
    with_t = normal_shock_condition(upstream_mach=3.0, upstream_temperature=300.0).rows[
        0
    ]
    assert old["velocity_ratio"] == with_t["velocity_ratio"]
    assert with_t["upstream_velocity"] == pytest.approx(3 * AIR.speed_of_sound(300.0))


@pytest.mark.parametrize(
    "model,temperature,missing",
    [("Keyes", 500.0, "downstream"), ("Blottner/Wilke", 500.0, "upstream")],
)
def test_transport_range_omits_only_affected_state(
    model: str, temperature: float, missing: str
) -> None:
    result = normal_shock_condition(
        upstream_mach=4.2,
        gas_model="HARMONIC_OSCILLATOR",
        upstream_temperature=temperature,
        upstream_pressure=25_000.0,
        viscosity_model=model,
        characteristic_length=0.5,
    )
    row = result.rows[0]
    assert row["status"] == "ok"
    assert row[f"{missing}_dynamic_viscosity"] is None
    assert row[f"{missing}_reynolds_number_per_length"] is None
    assert row[f"{missing}_reynolds_number"] is None
    present = "upstream" if missing == "downstream" else "downstream"
    assert number(row, f"{present}_reynolds_number_per_length") > 0
    assert number(row, "downstream_temperature") < 2000.0
    assert result.warnings
    assert "範囲外" in str(row["message"])


@pytest.mark.parametrize("kind", ["normal", "oblique", "cone"])
def test_sweeps_keep_transport_warnings_and_physical_failure_rows(kind: str) -> None:
    options: dict[str, Any] = dict(
        gas_model="HARMONIC_OSCILLATOR",
        upstream_temperature=500.0,
        upstream_pressure=25_000.0,
        viscosity_model="Blottner/Wilke",
    )
    if kind == "normal":
        result = normal_shock_sweep(start=3.0, stop=5.4, points=3, **options)
    elif kind == "oblique":
        result = oblique_shock_sweep(
            fixed_mach=3.0,
            fixed_deflection=0.0,
            branch=ShockBranch.WEAK,
            sweep_field="deflection",
            start=0.0,
            stop=float(np.deg2rad(60)),
            points=3,
            **options,
        )
    else:
        result = conical_shock_sweep(
            fixed_mach=3.0,
            fixed_cone_half_angle=0.0,
            sweep_field="cone_half_angle",
            start=0.0,
            stop=float(np.deg2rad(60)),
            points=3,
            **options,
        )
    assert result.rows[0]["status"] == "ok"
    assert result.rows[-1]["status"] in {"out_of_range", "no_attached_shock"}
    assert result.warnings
    assert result.rows[0]["upstream_dynamic_viscosity"] is None
    assert result.rows[-1].get("velocity_ratio") is None


@pytest.mark.parametrize("mach", [0.0, 2.0])
def test_isentropic_local_properties_without_pressure_and_transport_with_pressure(
    mach: float,
) -> None:
    common: dict[str, Any] = dict(
        input_value=mach,
        input_basis="mach",
        gas_model="NASA9",
        total_temperature=1000.0,
        with_heat_capacities=True,
        allow_extrapolation=False,
    )
    bare = isentropic_condition(**common).rows[0]
    result = isentropic_condition(
        total_pressure=50_000.0,
        viscosity_model="Keyes",
        characteristic_length=0.3,
        **common,
    )
    row = result.rows[0]
    t = number(bare, "static_temperature")
    assert number(bare, "heat_capacity_ratio") == pytest.approx(
        AIR_NASA9.heat_capacity_ratio(t), rel=2e-12
    )
    assert number(bare, "velocity") == pytest.approx(
        mach * float(AIR_NASA9.speed_of_sound(t)), rel=2e-12
    )
    assert bare["velocity"] == row["velocity"]
    assert bare.get("reynolds_number_per_length") is None
    expected = (
        number(row, "static_density")
        * number(row, "velocity")
        / float(AIR_KEYES_VISCOSITY.dynamic_viscosity(t))
    )
    assert row["reynolds_number_per_length"] == pytest.approx(expected, rel=2e-12)
    assert row["reynolds_number"] == pytest.approx(expected * 0.3, rel=2e-12)
    assert not result.warnings


@pytest.mark.parametrize("temperature", [79.0, 1845.0])
def test_transport_range_includes_both_endpoints(temperature: float) -> None:
    row = normal_shock_condition(
        upstream_mach=1.0,
        upstream_temperature=temperature,
        upstream_pressure=100_000.0,
        viscosity_model="Keyes",
    ).rows[0]
    assert number(row, "upstream_dynamic_viscosity") == number(
        row, "downstream_dynamic_viscosity"
    )


@pytest.mark.parametrize(
    "options,match",
    [
        ({"upstream_pressure": 0.0}, "pressure"),
        ({"upstream_pressure": float("nan")}, "pressure"),
        ({"characteristic_length": -1.0}, "characteristic_length"),
        ({"characteristic_length": 1.0}, "pressure is required"),
        ({"viscosity_model": "other"}, "viscosity_model"),
        ({"upstream_pressure": 1.0}, "upstream_temperature"),
    ],
)
def test_optional_input_validation(options: dict[str, object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        normal_shock_condition(upstream_mach=2.0, **options)  # type: ignore[arg-type]


def test_table_and_csv_optional_columns_and_inverse_length_units() -> None:
    row = normal_shock_condition(
        upstream_mach=3.0,
        upstream_temperature=500.0,
        upstream_pressure=25_000.0,
        with_heat_capacities=True,
        characteristic_length=0.2,
    ).rows[0]
    preferences = UnitPreferences(
        length="mm",
        inverse_length="1/mm",
        pressure="kPa",
        temperature="°C",
        speed="ft/s",
    )
    table = display_rows("normal_shock", (row,), preferences)
    assert table[0]["単位Re₁ [1/mm]"] == pytest.approx(
        number(row, "upstream_reynolds_number_per_length") / 1000
    )
    assert table[0]["静圧 p₁ [kPa]"] == 25.0
    assert table[0]["上流静温 T₁ [°C]"] == pytest.approx(226.85)
    assert "定圧比熱 cₚ₁ [J/(kg·K)]" in rows_to_csv(table)
    legacy = display_rows(
        "normal_shock", normal_shock_condition(upstream_mach=3.0).rows, preferences
    )
    assert "単位Re₁ [1/mm]" not in legacy[0]
    assert "定圧比熱 cₚ₁ [J/(kg·K)]" not in legacy[0]


def test_round_trip_temperature_at_upper_bound_keeps_valid_shock() -> None:
    from aerophysics.exceptions import ModelRangeError
    from aerophysics.shocks import normal_shock

    temperature = 404.40881763527057
    mach = 4.850018713943831
    core = normal_shock(mach, AIR_HARMONIC_OSCILLATOR, upstream_temperature=temperature)
    assert (2000.0 / temperature) * temperature > 2000.0
    assert temperature * float(core.static_temperature_ratio) == pytest.approx(
        2000.0, abs=2e-9
    )
    result = normal_shock_condition(
        upstream_mach=mach,
        gas_model="HARMONIC_OSCILLATOR",
        upstream_temperature=temperature,
        upstream_pressure=25_000.0,
    )
    row = result.rows[0]
    assert row["status"] == "ok"
    assert number(row, "downstream_temperature") <= 2000.0
    assert row["downstream_temperature"] == pytest.approx(2000.0, abs=2e-9)
    assert number(row, "downstream_heat_capacity_ratio") == pytest.approx(
        AIR_HARMONIC_OSCILLATOR.heat_capacity_ratio(2000.0), rel=2e-12
    )
    with pytest.raises(ModelRangeError):
        normal_shock_condition(
            upstream_mach=mach + 1e-9,
            gas_model="HARMONIC_OSCILLATOR",
            upstream_temperature=temperature,
            upstream_pressure=25_000.0,
        )


def test_temperature_restoration_preserves_one_ulp_and_extrapolation() -> None:
    from aerophysics.gui._flow_outputs import restore_static_temperature

    gas = AIR_HARMONIC_OSCILLATOR
    assert restore_static_temperature(float(np.nextafter(400.0, 0.0)), gas) == 400.0
    one_ulp = float(np.nextafter(2000.0, np.inf))
    two_ulps = float(np.nextafter(one_ulp, np.inf))
    assert restore_static_temperature(one_ulp, gas) == 2000.0
    assert restore_static_temperature(two_ulps, gas) == two_ulps
    assert restore_static_temperature(one_ulp, gas, allow_extrapolation=True) == one_ulp


def test_strong_branch_velocity_uses_full_flow_speed() -> None:
    row = oblique_shock_condition(
        upstream_mach=3.0,
        deflection_angle=float(np.deg2rad(10.0)),
        branch=ShockBranch.STRONG,
        gas_model="NASA9",
        upstream_temperature=500.0,
        upstream_pressure=25_000.0,
    ).rows[0]
    beta, theta = number(row, "shock_angle"), number(row, "deflection_angle")
    v1, v2 = number(row, "upstream_velocity"), number(row, "downstream_velocity")
    assert v1 * np.cos(beta) == pytest.approx(v2 * np.cos(beta - theta), rel=2e-10)
    assert row["velocity_ratio"] == pytest.approx(v2 / v1, rel=2e-12)
