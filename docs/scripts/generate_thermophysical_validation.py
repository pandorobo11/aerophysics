"""Generate thermophysical and real-gas verification artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from _verification_common import line_chart_svg, write_or_check

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "tests/reference_data/thermophysical"
TABLE_PATH = ROOT / "docs/_generated/thermophysical_validation.rst"
PROPERTY_PATH = ROOT / "docs/_static/thermophysical_properties.svg"
TRANSPORT_PATH = ROOT / "docs/_static/thermophysical_transport_accuracy.svg"
OFFICIAL_USSA = ROOT / "tests/reference_data/standard_atmosphere/official_1976.csv"


def _api() -> tuple[object, ...]:
    sys.path.insert(0, str(ROOT / "src"))
    from aerophysics.real_gas import (
        AIR_BEATTIE_BRIDGEMAN,
        AIR_HARMONIC_OSCILLATOR,
    )
    from aerophysics.thermochemistry import AIR_NASA7, AIR_NASA9
    from aerophysics.transport import (
        AIR_BLOTTNER_VISCOSITY,
        AIR_CONDUCTIVITY,
        AIR_KEYES_VISCOSITY,
        AIR_VISCOSITY,
    )

    return (
        AIR_NASA7,
        AIR_NASA9,
        AIR_HARMONIC_OSCILLATOR,
        AIR_BEATTIE_BRIDGEMAN,
        AIR_VISCOSITY,
        AIR_KEYES_VISCOSITY,
        AIR_BLOTTNER_VISCOSITY,
        AIR_CONDUCTIVITY,
    )


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _cantera_results(api: tuple[object, ...]) -> dict[str, tuple[float, float]]:
    gases = {"NASA7": api[0], "NASA9": api[1]}
    maxima = {
        name: (0.0, 0.0)
        for name in ("cp", "enthalpy", "entropy", "gamma", "speed of sound")
    }
    mapping = {
        "cp_J_kg_K": ("cp", "cp"),
        "enthalpy_J_kg": ("standard_enthalpy", "enthalpy"),
        "entropy_J_kg_K": ("entropy", "entropy"),
        "gamma": ("heat_capacity_ratio", "gamma"),
        "speed_of_sound_m_s": ("speed_of_sound", "speed of sound"),
    }
    for row in _rows(REFERENCE / "cantera-3.2.0.csv"):
        gas = gases[row["model"]]
        temperature = float(row["temperature_K"])
        pressure = float(row["pressure_Pa"])
        for column, (method_name, label) in mapping.items():
            method = getattr(gas, method_name)
            actual = (
                method(temperature, pressure)
                if method_name == "entropy"
                else method(temperature)
            )
            expected = float(row[column])
            relative = abs(float(actual) - expected) / max(abs(expected), 1.0)
            if relative > maxima[label][0]:
                maxima[label] = (relative, temperature)
    return maxima


def _invariants(api: tuple[object, ...]) -> list[tuple[str, str, float, bool]]:
    from aerophysics.isentropic import isentropic_state

    results: list[tuple[str, str, float, bool]] = []
    temperature = np.concatenate(
        (
            np.linspace(200.1, 999.0, 500),
            np.linspace(1001.0, 5999.9, 1000),
        )
    )
    for label, gas in (("NASA7", api[0]), ("NASA9", api[1])):
        cp = np.asarray(gas.cp(temperature))
        cv = np.asarray(gas.cv(temperature))
        gamma = np.asarray(gas.heat_capacity_ratio(temperature))
        gas_constant = gas.specific_gas_constant
        closure = max(
            float(np.max(np.abs((cp - cv) / gas_constant - 1.0))),
            float(np.max(np.abs(gamma / (cp / cv) - 1.0))),
        )
        step = 1.0e-2
        derivative = (
            np.asarray(gas.standard_enthalpy(temperature + step))
            - np.asarray(gas.standard_enthalpy(temperature - step))
        ) / (2.0 * step)
        derivative_error = float(np.max(np.abs(derivative / cp - 1.0)))
        boundary = 1000.0
        offset = 1.0e-6
        continuity = max(
            abs(
                float(gas.cp(boundary - offset)) / float(gas.cp(boundary + offset))
                - 1.0
            ),
            abs(
                float(gas.standard_enthalpy(boundary - offset))
                / float(gas.standard_enthalpy(boundary + offset))
                - 1.0
            ),
        )
        results.extend(
            [
                (
                    f"{label} heat-capacity identities",
                    "relative <= 1e-12",
                    closure,
                    closure <= 1.0e-12,
                ),
                (
                    f"{label} dh/dT = cp",
                    "relative <= 1e-7",
                    derivative_error,
                    derivative_error <= 1.0e-7,
                ),
                (
                    f"{label} 1000 K region continuity",
                    "relative <= 1e-7",
                    continuity,
                    continuity <= 1.0e-7,
                ),
            ]
        )
    harmonic = api[2]
    harmonic_temperature = np.linspace(400.0, 2000.0, 801)
    harmonic_closure = float(
        np.max(
            np.abs(
                np.asarray(harmonic.cp(harmonic_temperature))
                - np.asarray(harmonic.cv(harmonic_temperature))
                - harmonic.specific_gas_constant
            )
        )
        / harmonic.specific_gas_constant
    )
    results.append(
        (
            "Harmonic-oscillator cp-cv=R",
            "relative <= 1e-12",
            harmonic_closure,
            harmonic_closure <= 1.0e-12,
        )
    )
    real = api[3]
    minimum_dpdrho = float("inf")
    maximum_state_closure = 0.0
    for pressure in np.linspace(1.0e6, 10.0e6, 10):
        for value in np.linspace(400.0, 1200.0, 17):
            state = real.state(value, pressure)
            density = float(state.density)
            delta = density * 1.0e-5
            derivative = (
                float(real.pressure(value, density + delta))
                - float(real.pressure(value, density - delta))
            ) / (2.0 * delta)
            minimum_dpdrho = min(minimum_dpdrho, derivative)
            maximum_state_closure = max(
                maximum_state_closure,
                abs(float(state.pressure) / float(real.pressure(value, density)) - 1.0),
            )
    results.extend(
        [
            (
                "Beattie-Bridgeman stable density root",
                "minimum dp/drho > 0",
                minimum_dpdrho,
                minimum_dpdrho > 0.0,
            ),
            (
                "Beattie-Bridgeman pressure closure",
                "relative <= 1e-12",
                maximum_state_closure,
                maximum_state_closure <= 1.0e-12,
            ),
        ]
    )
    for label, gas, total_pressure, mach in (
        ("NASA7", api[0], 1.0e5, np.asarray([0.2, 1.0, 2.0])),
        ("NASA9", api[1], 1.0e5, np.asarray([0.2, 1.0, 2.0])),
        ("Harmonic oscillator", api[2], 1.0e5, np.asarray([0.2, 1.0, 2.0])),
        (
            "Beattie-Bridgeman",
            api[3],
            5.0e6,
            np.asarray([0.2, 1.0, 1.5]),
        ),
    ):
        total_temperature = 900.0 if label == "Beattie-Bridgeman" else 1500.0
        flow = isentropic_state(
            mach,
            gas,
            total_temperature=total_temperature,
            total_pressure=total_pressure,
        )
        static_temperature = np.asarray(flow.static_temperature)
        static_pressure = np.asarray(flow.static_pressure)
        velocity = np.asarray(flow.velocity)
        if label == "Beattie-Bridgeman":
            static_enthalpy = np.asarray(
                gas.enthalpy(static_temperature, static_pressure)
            )
            total_enthalpy = float(gas.enthalpy(total_temperature, total_pressure))
        else:
            static_enthalpy = np.asarray(gas.standard_enthalpy(static_temperature))
            total_enthalpy = float(gas.standard_enthalpy(total_temperature))
        static_entropy = np.asarray(gas.entropy(static_temperature, static_pressure))
        total_entropy = float(gas.entropy(total_temperature, total_pressure))
        enthalpy_closure = float(
            np.max(np.abs((static_enthalpy + 0.5 * velocity**2) / total_enthalpy - 1.0))
        )
        entropy_closure = float(
            np.max(
                np.abs(static_entropy - total_entropy) / max(abs(total_entropy), 1.0)
            )
        )
        results.extend(
            [
                (
                    f"{label} isentropic total enthalpy",
                    "relative <= 1e-12",
                    enthalpy_closure,
                    enthalpy_closure <= 1.0e-12,
                ),
                (
                    f"{label} isentropic entropy",
                    "relative <= 1e-12",
                    entropy_closure,
                    entropy_closure <= 1.0e-12,
                ),
            ]
        )
    return results


def _nist_assessment(
    api: tuple[object, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = _rows(REFERENCE / "nist-lemmon-jacobsen-2004.csv")
    temperature = np.asarray([float(row["temperature_K"]) for row in rows])
    viscosity = np.asarray([float(row["dynamic_viscosity_Pa_s"]) for row in rows])
    conductivity = np.asarray(
        [float(row["thermal_conductivity_W_m_K"]) for row in rows]
    )
    viscosity_difference = (
        np.asarray(api[4].dynamic_viscosity(temperature)) / viscosity - 1.0
    )
    conductivity_difference = (
        np.asarray(api[7].thermal_conductivity(temperature)) / conductivity - 1.0
    )
    return temperature, viscosity_difference, conductivity_difference


def _ussa_results(
    api: tuple[object, ...],
) -> list[tuple[str, str, float, float, bool]]:
    maximum_relative_difference = 0.0
    maximum_temperature = 0.0
    viscosity_passed = True
    for row in _rows(OFFICIAL_USSA):
        if not row["dynamic_viscosity_Pa_s"]:
            continue
        temperature = float(row["temperature_K"])
        expected = float(row["dynamic_viscosity_Pa_s"])
        printed_tolerance = float(row["dynamic_viscosity_abs_tolerance_Pa_s"])
        actual = float(api[4].dynamic_viscosity(temperature))
        difference = abs(actual - expected)
        tolerance = max(2.0 * printed_tolerance, 1.0e-4 * abs(expected))
        guard = 1.0e-14 * max(1.0, abs(expected))
        viscosity_passed &= difference <= tolerance + guard
        relative = difference / abs(expected)
        if relative > maximum_relative_difference:
            maximum_relative_difference = relative
            maximum_temperature = temperature

    metadata = json.loads(
        (REFERENCE / "ussa-1976-transport.json").read_text(encoding="utf-8")
    )
    conductivity_reference = metadata["thermal_conductivity"]
    conductivity_temperature = float(conductivity_reference["temperature_K"])
    conductivity_expected = float(conductivity_reference["value_W_m_K"])
    conductivity_tolerance = float(
        conductivity_reference["printed_abs_tolerance_W_m_K"]
    )
    conductivity_actual = float(api[7].thermal_conductivity(conductivity_temperature))
    conductivity_difference = abs(conductivity_actual - conductivity_expected)
    return [
        (
            "Sutherland / USSA Table III",
            "within 2 printed half-digits or 1e-4 relative",
            maximum_relative_difference,
            maximum_temperature,
            viscosity_passed,
        ),
        (
            "USSA conductivity / Equation (53) errata example",
            "absolute <= 5e-7 W/(m K)",
            conductivity_difference,
            conductivity_temperature,
            conductivity_difference <= conductivity_tolerance,
        ),
    ]


def _source_equation_results(
    api: tuple[object, ...],
) -> dict[str, tuple[float, float]]:
    models = {
        "Sutherland viscosity": api[4].dynamic_viscosity,
        "Keyes viscosity": api[5].dynamic_viscosity,
        "Blottner/Wilke viscosity": api[6].dynamic_viscosity,
        "USSA conductivity": api[7].thermal_conductivity,
    }
    maxima = {name: (0.0, 0.0) for name in models}
    for row in _rows(REFERENCE / "transport_source_equations.csv"):
        model = row["model"]
        temperature = float(row["temperature_K"])
        expected = float(row["value_SI"])
        actual = float(models[model](temperature))
        relative = abs(actual / expected - 1.0)
        if relative > maxima[model][0]:
            maxima[model] = (relative, temperature)
    return maxima


def _rst(
    cantera: dict[str, tuple[float, float]],
    ussa: list[tuple[str, str, float, float, bool]],
    transport: dict[str, tuple[float, float]],
    invariants: list[tuple[str, str, float, bool]],
    assessment: tuple[float, float],
) -> str:
    passed = (
        all(value[0] <= 2.0e-6 for value in cantera.values())
        and all(item[4] for item in ussa)
        and all(value[0] <= 1.0e-12 for value in transport.values())
        and all(item[3] for item in invariants)
    )
    status = "Verified" if passed else "Needs revision"
    lines = [
        f"**Overall status: {status}.**",
        "",
        "Cantera 3.2.0 comparison",
        "~~~~~~~~~~~~~~~~~~~~~~~~",
        "",
        ".. list-table:: Maximum relative differences from the pinned snapshot",
        "   :header-rows: 1",
        "",
        "   * - Quantity",
        "     - Maximum relative difference",
        "     - Temperature [K]",
        "     - Result",
    ]
    for label, (difference, temperature) in cantera.items():
        lines.extend(
            [
                f"   * - {label}",
                f"     - {difference:.4g}",
                f"     - {temperature:g}",
                f"     - {'Pass' if difference <= 2.0e-6 else 'Fail'}",
            ]
        )
    lines.extend(
        [
            "",
            "Primary transport references",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            ".. list-table:: Comparisons with U.S. Standard Atmosphere 1976",
            "   :header-rows: 1",
            "",
            "   * - Reference",
            "     - Criterion",
            "     - Maximum difference",
            "     - Temperature [K]",
            "     - Result",
        ]
    )
    for label, criterion, difference, temperature, item_passed in ussa:
        lines.extend(
            [
                f"   * - {label}",
                f"     - {criterion}",
                f"     - {difference:.4g}",
                f"     - {temperature:g}",
                f"     - {'Pass' if item_passed else 'Fail'}",
            ]
        )
    lines.extend(
        [
            "",
            "Published transport equations",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            ".. list-table:: Maximum relative differences from direct "
            "source-equation reproductions",
            "   :header-rows: 1",
            "",
            "   * - Model",
            "     - Maximum relative difference",
            "     - Temperature [K]",
            "     - Result",
        ]
    )
    for label, (difference, temperature) in transport.items():
        lines.extend(
            [
                f"   * - {label}",
                f"     - {difference:.4g}",
                f"     - {temperature:g}",
                f"     - {'Pass' if difference <= 1.0e-12 else 'Fail'}",
            ]
        )
    lines.extend(
        [
            "",
            "NIST physical-accuracy assessment",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            (
                "The Lemmon--Jacobsen evaluated dilute-air correlation is an "
                "independent physical-accuracy reference, not an acceptance test "
                "for the intentionally simpler Sutherland and USSA correlations. "
                "Across 250--1500 K at zero density, the largest absolute relative "
                f"differences are ``{assessment[0]:.3%}`` for viscosity and "
                f"``{assessment[1]:.3%}`` for conductivity. The NIST source "
                "estimates dilute-gas uncertainties of 1% and 2%, respectively, "
                "over this range."
            ),
            "",
            "Thermodynamic invariants",
            "~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            ".. list-table:: Dense-grid checks",
            "   :header-rows: 1",
            "",
            "   * - Check",
            "     - Criterion",
            "     - Maximum or minimum",
            "     - Result",
        ]
    )
    for label, criterion, value, item_passed in invariants:
        lines.extend(
            [
                f"   * - {label}",
                f"     - {criterion}",
                f"     - {value:.4g}",
                f"     - {'Pass' if item_passed else 'Fail'}",
            ]
        )
    return "\n".join(lines)


def generate(*, check: bool) -> bool:
    api = _api()
    cantera = _cantera_results(api)
    ussa_results = _ussa_results(api)
    transport_results = _source_equation_results(api)
    invariants = _invariants(api)
    temperature, viscosity_difference, conductivity_difference = _nist_assessment(api)
    assessment = (
        float(np.max(np.abs(viscosity_difference))),
        float(np.max(np.abs(conductivity_difference))),
    )
    profile_temperature = np.linspace(200.0, 6000.0, 500)
    cp7 = np.asarray(api[0].cp(profile_temperature))
    cp9 = np.asarray(api[1].cp(profile_temperature))
    properties = line_chart_svg(
        title="Frozen dry-air heat capacity",
        description=(
            "NASA seven- and nine-coefficient frozen-air constant-pressure heat "
            "capacities from 200 to 6000 kelvin."
        ),
        x_label="Temperature [K]",
        y_label="cp [J/(kg K)]",
        x=profile_temperature,
        series=(("NASA7", "#0072B2", cp7), ("NASA9", "#D55E00", cp9)),
    )
    transport = line_chart_svg(
        title="Transport-model differences from NIST evaluated dilute air",
        description=(
            "Non-gating physical-accuracy differences for the Sutherland viscosity "
            "and USSA conductivity correlations relative to Lemmon and Jacobsen's "
            "evaluated dilute-air reference."
        ),
        x_label="Temperature [K]",
        y_label="aerophysics / NIST reference - 1 [-]",
        x=temperature,
        series=(
            ("Viscosity", "#009E73", viscosity_difference),
            ("Conductivity", "#CC79A7", conductivity_difference),
        ),
    )
    current = True
    current &= write_or_check(
        TABLE_PATH,
        _rst(cantera, ussa_results, transport_results, invariants, assessment),
        check=check,
    )
    current &= write_or_check(PROPERTY_PATH, properties, check=check)
    current &= write_or_check(TRANSPORT_PATH, transport, check=check)
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not generate(check=args.check):
        raise SystemExit("thermophysical verification artifacts are stale")


if __name__ == "__main__":
    main()
