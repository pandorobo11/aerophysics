"""Generate boundary-layer and protrusion verification artifacts."""

from __future__ import annotations

import argparse
import csv
import sys
import warnings
from pathlib import Path

import numpy as np
from _verification_common import format_bounded_error, line_chart_svg, write_or_check
from scipy.integrate import cumulative_trapezoid, solve_bvp

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "tests/reference_data/viscous_flow"
TABLE_PATH = ROOT / "docs/_generated/viscous_flow_validation.rst"
FRICTION_PATH = ROOT / "docs/_static/viscous_skin_friction.svg"
PROTRUSION_PATH = ROOT / "docs/_static/protrusion_shielding.svg"


def _api() -> tuple[object, ...]:
    sys.path.insert(0, str(ROOT / "src"))
    from aerophysics.boundary_layer import (
        BoundaryLayerRegime,
        TurbulentModel,
        flat_plate_boundary_layer,
    )
    from aerophysics.boundary_layer_profile import (
        CompressibleVelocityTransformation,
        transform_compressible_velocity_profile,
    )
    from aerophysics.gas import AIR
    from aerophysics.protrusion import protrusion_drag
    from aerophysics.transport import AIR_VISCOSITY

    return (
        BoundaryLayerRegime,
        TurbulentModel,
        CompressibleVelocityTransformation,
        flat_plate_boundary_layer,
        transform_compressible_velocity_profile,
        protrusion_drag,
        AIR_VISCOSITY,
        AIR,
    )


def _rows() -> list[dict[str, str]]:
    with (REFERENCE / "flat_plate_source_equations.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        return list(csv.DictReader(stream))


def _source_comparison(api: tuple[object, ...]) -> dict[str, tuple[float, float]]:
    regime_enum, model_enum, _, flat_plate, _, _, _, _ = api
    maxima = {
        name: (0.0, 0.0)
        for name in (
            "Blasius delta99",
            "Blasius displacement thickness",
            "Blasius momentum thickness",
            "Blasius local Cf",
            "Blasius average Cf",
            "one-fifth-power thickness delta99",
            "Schlichting average Cf",
        )
    }
    for row in _rows():
        reynolds = float(row["reynolds_number"])
        laminar = flat_plate(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=regime_enum.LAMINAR,
        )
        schlichting = flat_plate(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=regime_enum.TURBULENT,
            turbulent_model=model_enum.SCHLICHTING,
        )
        comparisons = {
            "Blasius delta99": (
                laminar.boundary_layer_thickness,
                "blasius_delta_over_x",
            ),
            "Blasius displacement thickness": (
                laminar.displacement_thickness,
                "blasius_displacement_over_x",
            ),
            "Blasius momentum thickness": (
                laminar.momentum_thickness,
                "blasius_momentum_over_x",
            ),
            "Blasius local Cf": (
                laminar.local_skin_friction_coefficient,
                "blasius_local_cf",
            ),
            "Blasius average Cf": (
                laminar.average_skin_friction_coefficient,
                "blasius_average_cf",
            ),
            "one-fifth-power thickness delta99": (
                schlichting.boundary_layer_thickness,
                "power_delta_over_x",
            ),
            "Schlichting average Cf": (
                schlichting.average_skin_friction_coefficient,
                "schlichting_average_cf",
            ),
        }
        for label, (actual, column) in comparisons.items():
            expected = float(row[column])
            relative = abs(float(actual) / expected - 1.0)
            if relative > maxima[label][0]:
                maxima[label] = (relative, reynolds)
    return maxima


def _chart_comparison(api: tuple[object, ...]) -> tuple[float, float]:
    regime_enum, model_enum, _, flat_plate, _, _, viscosity_model, gas = api
    maximum = (0.0, 0.0)
    with (REFERENCE / "van_driest_ii_chart.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        for row in csv.DictReader(stream):
            temperature = float(row["edge_temperature_K"])
            mach = float(row["mach"])
            reynolds = float(row["reynolds_number"])
            edge_viscosity = float(viscosity_model.dynamic_viscosity(temperature))
            recovery = temperature * (
                1.0 + np.cbrt(0.72) * (gas.heat_capacity_ratio - 1.0) * mach**2 / 2.0
            )
            result = flat_plate(
                1.0,
                reynolds * edge_viscosity,
                1.0,
                edge_viscosity,
                regime=regime_enum.TURBULENT,
                compressible=True,
                turbulent_model=model_enum.VAN_DRIEST_II,
                mach=mach,
                edge_temperature=temperature,
                wall_temperature=(
                    float(row["wall_to_adiabatic_temperature_ratio"]) * recovery
                ),
            )
            ratio = abs(
                float(result.local_skin_friction_coefficient) - float(row["local_cf"])
            ) / float(row["absolute_tolerance"])
            if ratio > maximum[0]:
                maximum = (ratio, mach)
    return maximum


def _invariants(api: tuple[object, ...]) -> list[tuple[str, str, float, bool]]:
    (
        regime_enum,
        model_enum,
        transformation_enum,
        flat_plate,
        transform_profile,
        protrusion_drag,
        viscosity_model,
        _,
    ) = api
    reynolds = np.geomspace(1.0e5, 1.0e9, 300)
    result = flat_plate(
        1.0,
        reynolds,
        1.0,
        1.0,
        regime=regime_enum.TURBULENT,
    )
    average = np.asarray(result.average_skin_friction_coefficient)
    expected_local = average * (1.0 - 2.58 / np.log(reynolds))
    derivative_error = float(
        np.max(
            np.abs(
                np.asarray(result.local_skin_friction_coefficient) / expected_local
                - 1.0
            )
        )
    )
    edge_viscosity = float(viscosity_model.dynamic_viscosity(300.0))
    incompressible = flat_plate(
        1.0,
        100.0,
        1.0,
        edge_viscosity,
        regime=regime_enum.TURBULENT,
        compressible=True,
        turbulent_model=model_enum.VAN_DRIEST_II,
        mach=0.0,
        edge_temperature=300.0,
    )
    low_mach = flat_plate(
        1.0,
        100.0,
        1.0,
        edge_viscosity,
        regime=regime_enum.TURBULENT,
        compressible=True,
        turbulent_model=model_enum.VAN_DRIEST_II,
        mach=1.0e-6,
        edge_temperature=300.0,
    )
    low_mach_error = abs(
        float(low_mach.average_skin_friction_coefficient)
        / float(incompressible.average_skin_friction_coefficient)
        - 1.0
    )
    height = np.linspace(0.0, 0.01, 1001)
    velocity = np.linspace(0.0, 100.0, 1001)
    transformed = transform_profile(
        height,
        velocity,
        np.full_like(height, 1.2),
        np.full_like(height, 1.8e-5),
        12.0,
        transformation=transformation_enum.VAN_DRIEST,
    )
    friction_velocity = np.sqrt(12.0 / 1.2)
    profile_error = max(
        float(
            np.max(
                np.abs(
                    np.asarray(transformed.transformed_wall_coordinate)
                    - height * 1.2 * friction_velocity / 1.8e-5
                )
            )
        ),
        float(
            np.max(
                np.abs(
                    np.asarray(transformed.transformed_velocity_plus)
                    - velocity / friction_velocity
                )
            )
        ),
    )
    height_ratio = np.asarray([0.01, 0.03, 0.1, 0.3, 0.8])
    shielding_error = 0.0
    for ratio in height_ratio:
        drag = protrusion_drag(
            1.0,
            float(ratio),
            1.0,
            1.0,
            1.0,
            1.0,
            integration_points=32_769,
        )
        expected = 7.0 / 9.0 * ratio ** (2.0 / 7.0)
        shielding_error = max(
            shielding_error, abs(float(drag.shielding_factor) / expected - 1.0)
        )
    constant = protrusion_drag(
        1.2,
        0.02,
        0.01,
        100.0,
        1.2,
        0.1,
        profile_height=[0.0, 0.02],
        profile_velocity=[100.0, 100.0],
        profile_density=[1.2, 1.2],
    )
    expected_force = 1.2 * 0.5 * 1.2 * 100.0**2 * 0.02 * 0.01
    constant_error = abs(float(constant.direct_drag) / expected_force - 1.0)
    coarse = protrusion_drag(1.0, 0.3, 1.0, 1.0, 1.0, 1.0, integration_points=4097)
    fine = protrusion_drag(1.0, 0.3, 1.0, 1.0, 1.0, 1.0, integration_points=32_769)
    grid_error = abs(
        float(coarse.shielding_factor) / float(fine.shielding_factor) - 1.0
    )
    vanishing = protrusion_drag(
        1.0, 1.0e-12, 1.0, 1.0, 1.0, 1.0, integration_points=4097
    )
    outside = protrusion_drag(1.0, 1.0e4, 1.0, 1.0, 1.0, 1.0, integration_points=4097)
    return [
        (
            "Schlichting local/average derivative",
            "relative <= 1e-12",
            derivative_error,
            derivative_error <= 1.0e-12,
        ),
        (
            "Van Driest II low-Mach limit",
            "relative <= 1e-10",
            low_mach_error,
            low_mach_error <= 1.0e-10,
        ),
        (
            "Van Driest constant-property transform",
            "absolute <= 1e-10",
            profile_error,
            profile_error <= 1.0e-10,
        ),
        (
            "One-seventh-power protrusion integral",
            "relative <= 1e-6",
            shielding_error,
            shielding_error <= 1.0e-6,
        ),
        (
            "Constant-profile protrusion integral",
            "relative <= 1e-12",
            constant_error,
            constant_error <= 1.0e-12,
        ),
        (
            "Protrusion integration grid convergence",
            "relative <= 1e-5",
            grid_error,
            grid_error <= 1.0e-5,
        ),
        (
            "Protrusion zero-height limit",
            "shielding < 1e-3",
            float(vanishing.shielding_factor),
            float(vanishing.shielding_factor) < 1.0e-3,
        ),
        (
            "Protrusion outside-layer limit",
            "absolute distance from 1 <= 1e-4",
            abs(float(outside.shielding_factor) - 1.0),
            abs(float(outside.shielding_factor) - 1.0) <= 1.0e-4,
        ),
    ]


def _rst(
    source: dict[str, tuple[float, float]],
    chart: tuple[float, float],
    invariants: list[tuple[str, str, float, bool]],
) -> str:
    passed = (
        all(value[0] <= 1.0e-12 for value in source.values())
        and chart[0] <= 1.0
        and all(item[3] for item in invariants)
    )
    status = "Verified" if passed else "Needs revision"
    lines = [
        f"**Overall status: {status}.**",
        "",
        "Published-correlation comparison",
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
        "",
        (
            ".. list-table:: Maximum relative differences from independently "
            "evaluated source equations"
        ),
        "   :header-rows: 1",
        "",
        "   * - Quantity",
        "     - Maximum relative difference",
        "     - Reynolds number",
        "     - Result",
    ]
    for label, (difference, reynolds) in source.items():
        lines.extend(
            [
                f"   * - {label}",
                f"     - {difference:.4g}",
                f"     - {reynolds:.4g}",
                f"     - {'Pass' if difference <= 1.0e-12 else 'Fail'}",
            ]
        )
    lines.extend(
        [
            "",
            "NASA TN D-6945 chart comparison",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            (
                "For representative Figure 3(e) points, the largest absolute "
                f"difference/tolerance ratio is ``{chart[0]:.4g}`` at Mach "
                f"``{chart[1]:g}`` (``{'Pass' if chart[0] <= 1.0 else 'Fail'}``)."
            ),
        ]
    )
    lines.extend(
        [
            "",
            "Physical and numerical invariants",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            ".. list-table:: Profile, limit, and integration checks",
            "   :header-rows: 1",
            "",
            "   * - Check",
            "     - Criterion",
            "     - Maximum difference",
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


def _laminar_similarity_accuracy(api: tuple[object, ...]) -> str:
    regime_enum, _, _, flat_plate, _, _, viscosity_model, _ = api
    with (REFERENCE / "laminar_compressible_similarity.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]
    viscosity = float(viscosity_model.dynamic_viscosity(240.0))
    errors = []
    for row in rows:
        result = flat_plate(
            1.0,
            1e8 * viscosity,
            1.0,
            viscosity,
            regime=regime_enum.LAMINAR,
            compressible=True,
            mach=row["mach"],
            edge_temperature=240.0,
            wall_temperature=(
                None
                if row["wall_to_adiabatic_temperature_ratio"] == 1
                else row["wall_temperature_K"]
            ),
        )
        errors.append(
            [
                100.0 * (float(value) * 1e4 / row[key] - 1.0)
                for value, key in (
                    (result.local_skin_friction_coefficient, "cf_sqrt_re"),
                    (result.momentum_thickness, "momentum_sqrt_re_over_x"),
                    (result.boundary_layer_thickness, "delta99_sqrt_re_over_x"),
                    (result.displacement_thickness, "displacement_sqrt_re_over_x"),
                )
            ]
        )
    differences = np.asarray(errors)
    displacement = np.array([row["displacement_sqrt_re_over_x"] for row in rows])
    relative = np.abs(displacement) >= 0.05
    absolute_errors = differences[:, 3] / 100 * displacement
    lines = [
        ".. list-table:: Library versus independent laminar similarity",
        "   :header-rows: 1",
        "",
        "   * - Quantity / subset",
        "     - Points",
        "     - Error unit",
        "     - Maximum absolute",
        "     - RMS",
    ]
    for label, values, unit in (
        ("Skin friction", differences[:, 0], "%"),
        ("Momentum thickness", differences[:, 1], "%"),
        ("99-percent thickness", differences[:, 2], "%"),
        ("Displacement, coefficient magnitude >= 0.05", differences[relative, 3], "%"),
        (
            "Displacement, coefficient magnitude < 0.05",
            absolute_errors[~relative],
            "coefficient",
        ),
    ):
        limit = 1e-9 if unit == "coefficient" else 1e-7
        maximum = format_bounded_error(float(np.max(np.abs(values))), limit, digits=3)
        rms = format_bounded_error(float(np.sqrt(np.mean(values**2))), limit, digits=3)
        lines.extend(
            [
                f"   * - {label}",
                f"     - {len(values)}",
                f"     - {unit}",
                f"     - {maximum}",
                f"     - {rms}",
            ]
        )
    return "\n".join(lines)


def _turbulent_sa_comparison(api: tuple[object, ...]) -> str:
    regime_enum, model_enum, _, flat_plate, _, _, viscosity_model, gas = api
    with (REFERENCE / "turbulent_sa_boundary_layer.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = list(csv.DictReader(stream))
    lines = [
        ".. list-table:: ``VAN_DRIEST_II_WALL_TEMPERATURE`` minus the in-repository "
        "SA boundary-layer reference, divided by the reference",
        "   :header-rows: 1",
        "",
        "   * - Case",
        "     - :math:`T_w/T_r`",
        "     - :math:`Re_{x,e}`",
        "     - :math:`C_f`",
        "     - :math:`\\delta_{99}`",
        "     - :math:`\\delta^*`",
        "     - :math:`\\theta`",
        "     - :math:`q_w`",
    ]
    for row in rows:
        temperature = float(row["edge_temperature_K"])
        reynolds = float(row["reynolds_number"])
        viscosity = float(viscosity_model.dynamic_viscosity(temperature))
        # The profile closure uses the edge velocity, so it must match Me.
        velocity = float(row["mach"]) * np.sqrt(
            gas.heat_capacity_ratio * gas.specific_gas_constant * temperature
        )
        result = flat_plate(
            1.0,
            velocity,
            reynolds * viscosity / velocity,
            viscosity,
            regime=regime_enum.TURBULENT,
            compressible=True,
            turbulent_model=model_enum.VAN_DRIEST_II_WALL_TEMPERATURE,
            mach=float(row["mach"]),
            edge_temperature=temperature,
            wall_temperature=float(row["wall_temperature_K"]),
        )
        differences = [
            100.0 * (float(value) / float(row[key]) - 1.0)
            for value, key in (
                (result.local_skin_friction_coefficient, "cf"),
                (result.boundary_layer_thickness, "delta99_over_x"),
                (result.displacement_thickness, "displacement99_over_x"),
                (result.momentum_thickness, "momentum_over_x"),
            )
        ]
        ratio = float(row["wall_to_recovery_temperature_ratio"])
        # An adiabatic correlation wall has no flux to compare.
        heat_flux = "n/a"
        if ratio != 1.0:
            coefficient = float(result.wall_heat_flux) / (
                reynolds * viscosity * gas.cp * temperature
            )
            expected = float(row["wall_heat_flux_coefficient"])
            heat_flux = f"{100.0 * (coefficient / expected - 1.0):+.1f}%"
        lines.extend(
            [
                f"   * - {row['case']}",
                f"     - {ratio:.2f}",
                f"     - {reynolds:.0e}",
                *(f"     - {value:+.1f}%" for value in differences),
                f"     - {heat_flux}",
            ]
        )
    return "\n".join(lines)


def _similarity_reference(beta: float, wall_ratio: float) -> tuple[float, float]:
    """Solve the similar equations with ``C = Pr = 1`` and no dissipation."""

    def rhs(_: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.vstack(
            (y[1], y[2], -y[0] * y[2] - beta * (y[3] - y[1] ** 2), y[4], -y[0] * y[4])
        )

    def boundary(wall: np.ndarray, edge: np.ndarray) -> np.ndarray:
        return np.array(
            (wall[0], wall[1], edge[1] - 1.0, wall[3] - wall_ratio, edge[3] - 1.0)
        )

    x = np.linspace(0.0, 12.0, 400)
    decay = np.exp(-x)
    guess = np.vstack(
        (
            x - 1.0 + decay,
            1.0 - decay,
            decay,
            1.0 + (wall_ratio - 1.0) * decay,
            (1.0 - wall_ratio) * decay,
        )
    )
    solution = solve_bvp(rhs, boundary, x, guess, tol=1e-10, max_nodes=200000)
    return float(solution.y[2, 0]), float(solution.y[4, 0])


class _LinearViscosity:
    """Viscosity proportional to temperature, which makes ``C = 1``."""

    def dynamic_viscosity(self, temperature: object) -> object:
        value = np.asarray(temperature, dtype=np.float64) * 1e-7
        return float(value) if np.ndim(temperature) == 0 else value


def _howarth_separation(mach0: float, marching: tuple[object, ...]) -> float:
    from aerophysics.gas import AIR

    edge_type, regime, _, separation_error = marching
    gamma, gas_constant = AIR.heat_capacity_ratio, AIR.specific_gas_constant
    transformed = np.linspace(0.0, 1.0, 161)
    mach = mach0 * (1.0 - transformed / 8.0)
    temperature = 300.0 / (1.0 + (gamma - 1.0) / 2.0 * mach**2)
    pressure = 1e5 * (temperature / 300.0) ** (gamma / (gamma - 1.0))
    density = pressure / (gas_constant * temperature)
    velocity = mach * np.sqrt(gamma * gas_constant * temperature)
    model = _LinearViscosity()
    distance = cumulative_trapezoid(
        mach / (density * model.dynamic_viscosity(temperature) * velocity),
        transformed,
        initial=0,
    )
    distance /= distance[-1]
    edge = edge_type.from_pressure(
        distance,
        pressure,
        total_pressure=1e5,
        total_temperature=300.0,
        viscosity_model=model,
    )
    try:
        marching[2](edge, regime=regime.LAMINAR, prandtl_number=1.0)
    except separation_error as caught:
        position = float(np.interp(caught.separation_distance, distance, transformed))
        return position / 8.0
    raise RuntimeError("the Howarth case did not separate")


def _integral_balances(marching: tuple[object, ...]) -> tuple[float, float]:
    from aerophysics.gas import AIR

    edge_type, regime, solve, _ = marching
    gamma, gas_constant = AIR.heat_capacity_ratio, AIR.specific_gas_constant
    cp = AIR.cp
    distance = np.linspace(0.0, 2.0, 201)
    mach = 2.0 + (1.0 + np.tanh((distance - 1.0) / 0.1)) / 2.0
    pressure = 2e5 * (1.0 + (gamma - 1.0) / 2.0 * mach**2) ** (-gamma / (gamma - 1.0))
    edge = edge_type.from_pressure(
        distance, pressure, total_pressure=2e5, total_temperature=450.0
    )
    result = solve(edge, regime=regime.LAMINAR, wall_temperature=300.0)
    density = edge.pressure[1:] / (gas_constant * edge.temperature[1:])
    velocity = edge.velocity[1:]
    stations = result.distance
    momentum = (
        np.gradient(density * velocity**2 * result.momentum_thickness, stations)
        + density
        * velocity
        * np.gradient(velocity, stations)
        * result.displacement_thickness
    ) / result.wall_shear_stress
    deficit = np.array(
        [
            np.trapezoid(
                density[i]
                / result.temperature_ratio[i]
                * velocity[i]
                * result.velocity_ratio[i]
                * (
                    cp * edge.temperature[i + 1] * (result.temperature_ratio[i] - 1.0)
                    + velocity[i] ** 2 / 2.0 * (result.velocity_ratio[i] ** 2 - 1.0)
                ),
                result.wall_distance[i],
            )
            for i in range(stations.size)
        ]
    )
    energy = -np.gradient(deficit, stations) / result.wall_heat_flux
    mask = stations > 0.05
    return (
        float(np.median(np.abs(momentum[mask] - 1.0))),
        float(np.median(np.abs(energy[mask] - 1.0))),
    )


def _marching_laminar_verification(api: tuple[object, ...]) -> str:
    from aerophysics._boundary_layer_march import StationEdge, eta_grid, march
    from aerophysics.exceptions import BoundaryLayerSeparationError
    from aerophysics.marching_boundary_layer import (
        BoundaryLayerEdge,
        _EdgeInterpolant,
        _result,
        marching_boundary_layer,
    )

    regime, _, _, flat_plate, _, _, viscosity_model, gas = api
    gamma, gas_constant = gas.heat_capacity_ratio, gas.specific_gas_constant
    marching = (
        BoundaryLayerEdge,
        regime,
        marching_boundary_layer,
        BoundaryLayerSeparationError,
    )
    eta = eta_grid(1e-3, 1.015, 0.05, 16.0)

    def wedge(beta: float, wall_ratio: float) -> tuple[float, float]:
        solution = march(
            lambda xi: StationEdge(xi, 1.0, beta, 0.0, 300.0, wall_ratio),
            np.array([1.0]),
            eta=eta,
            prandtl_number=1.0,
            viscosity_model=_LinearViscosity(),
        )
        return (
            float(solution.wall_velocity_slope[-1]),
            float(solution.wall_temperature_slope[-1]),
        )

    rows: list[tuple[str, str, float, float]] = []
    for beta in (-0.15, 0.0, 0.5):
        computed, _ = wedge(beta, 1.0)
        reference, _ = _similarity_reference(beta, 1.0)
        rows.append((f"Falkner--Skan, beta = {beta:g}", "f''(0)", computed, reference))
    for beta, wall_ratio in ((-0.1, 0.2), (0.5, 2.0)):
        shear, heat = wedge(beta, wall_ratio)
        reference_shear, reference_heat = _similarity_reference(beta, wall_ratio)
        label = f"Similar, beta = {beta:g}, g_w = {wall_ratio:g}"
        rows.append((label, "f''(0)", shear, reference_shear))
        rows.append((label, "g'(0)", heat, reference_heat))
    for beta, wall_enthalpy, published_shear, published_heat in (
        (1.0, 0.0, 1.2326, None),
        (2.0, 0.0, 1.6870, None),
        (0.5, 1.0, 1.2351, -0.5725),
        (2.0, 1.0, 2.4878, -0.6613),
    ):
        shear, heat = wedge(beta, 1.0 + wall_enthalpy)
        label = f"NACA 1293 Table II, beta = {beta:g}, S_w = {wall_enthalpy:g}"
        rows.append((label, "f''(0)", shear, published_shear))
        if published_heat is not None:
            rows.append((label, "g'(0)", heat, published_heat))
    incompressible = _howarth_separation(0.01, marching)
    rows.append(("Howarth flow, M = 0.01", "x_sep/L", incompressible, 0.9589 / 8.0))
    rows.append(
        (
            "Howarth image, M0 = 3",
            "x_sep/L",
            _howarth_separation(3.0, marching),
            incompressible,
        )
    )
    for mach, temperature, wall in ((2.0, 220.0, None), (10.0, 60.0, 300.0)):
        pressure = 2e4
        ratio = 1.0 + (gamma - 1.0) / 2.0 * mach**2
        edge = BoundaryLayerEdge.from_pressure(
            [0.0, 1.0],
            [pressure, pressure],
            total_pressure=pressure * ratio ** (gamma / (gamma - 1.0)),
            total_temperature=temperature * ratio,
        )
        density = pressure / (gas_constant * temperature)
        velocity = mach * np.sqrt(gamma * gas_constant * temperature)
        viscosity = float(viscosity_model.dynamic_viscosity(temperature))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = marching_boundary_layer(
                edge, regime=regime.LAMINAR, wall_temperature=wall
            )
            plate = flat_plate(
                1.0,
                velocity,
                density,
                viscosity,
                regime=regime.LAMINAR,
                compressible=True,
                mach=mach,
                edge_temperature=temperature,
                wall_temperature=wall,
            )
        label = f"Plate, M = {mach:g}, " + (
            "adiabatic" if wall is None else f"T_w = {wall:g} K"
        )
        rows.append(
            (
                label,
                "C_f",
                float(result.local_skin_friction_coefficient[-1]),
                float(plate.local_skin_friction_coefficient),
            )
        )
        rows.append(
            (
                label,
                "theta",
                float(result.momentum_thickness[-1]),
                float(plate.momentum_thickness),
            )
        )
    wall_change = np.array([0.0, 0.5, 0.51, 1.0])
    wall_edge = BoundaryLayerEdge.from_pressure(
        wall_change,
        np.full(4, 2e4),
        total_pressure=1.6e5,
        total_temperature=500.0,
    )
    wall = np.array([300.0, 300.0, 400.0, 400.0])
    wall_result = marching_boundary_layer(
        wall_edge, regime=regime.LAMINAR, wall_temperature=wall
    )
    interpolant = _EdgeInterpolant(wall_edge, wall)
    stations = np.unique(np.concatenate((wall_change[1:], np.linspace(0.5, 0.51, 161))))
    refined = _result(
        wall_edge,
        interpolant,
        march(
            interpolant.state,
            stations,
            eta=eta,
            prandtl_number=0.72,
            viscosity_model=wall_edge.viscosity_model,
        ),
        0.72,
    )
    rows.append(
        (
            "Wall rises 100 K over 1 cm, s = 0.51 m",
            "q_w",
            float(wall_result.wall_heat_flux[1]),
            float(refined.wall_heat_flux[np.argmin(np.abs(stations - 0.51))]),
        )
    )
    momentum, energy = _integral_balances(marching)
    lines = [
        ".. list-table:: Marching laminar boundary-layer verification",
        "   :header-rows: 1",
        "   :widths: 40 15 15 15 15",
        "",
        "   * - Case",
        "     - Quantity",
        "     - Computed",
        "     - Reference",
        "     - Difference [%]",
    ]
    for label, quantity, computed, reference in rows:
        lines.extend(
            [
                f"   * - {label}",
                f"     - ``{quantity}``",
                f"     - {computed:.4g}",
                f"     - {reference:.4g}",
                f"     - {100.0 * (computed / reference - 1.0):+.2f}",
            ]
        )
    for label, value in (("Momentum", momentum), ("Total enthalpy", energy)):
        lines.extend(
            [
                f"   * - {label} integral, M = 2 to 3, T_w = 300 K",
                "     - median residual",
                "     - --",
                "     - --",
                f"     - {100.0 * value:.2f}",
            ]
        )
    return "\n".join(lines)


def generate(*, check: bool) -> bool:
    api = _api()
    source = _source_comparison(api)
    chart = _chart_comparison(api)
    invariants = _invariants(api)
    reynolds = np.geomspace(1.0e5, 1.0e9, 400)
    laminar = 0.664 / np.sqrt(reynolds)
    schlichting_average = 0.455 / np.log10(reynolds) ** 2.58
    schlichting = schlichting_average * (1.0 - 2.58 / np.log(reynolds))
    friction = line_chart_svg(
        title="Smooth flat-plate local skin friction",
        description=(
            "Blasius and Schlichting local skin-friction "
            "correlations over Reynolds number."
        ),
        x_label="Reynolds number [-]",
        y_label="Local Cf [-]",
        x=reynolds,
        series=(
            ("Blasius", "#0072B2", laminar),
            ("Schlichting", "#009E73", schlichting),
        ),
        log_x=True,
        log_y=True,
    )
    ratio = np.geomspace(1.0e-3, 1.0, 250)
    shielding = 7.0 / 9.0 * ratio ** (2.0 / 7.0)
    protrusion = line_chart_svg(
        title="Boundary-layer shielding of a protrusion",
        description=(
            "Closed-form effective dynamic-pressure factor for a constant-width "
            "protrusion in a one-seventh-power velocity profile."
        ),
        x_label="Protrusion height / boundary-layer thickness [-]",
        y_label="Shielding factor [-]",
        x=ratio,
        series=(("Analytic 7/9 (h/delta)^(2/7)", "#CC79A7", shielding),),
        log_x=True,
    )
    current = True
    current &= write_or_check(
        ROOT / "docs/_generated/laminar_similarity_accuracy.rst",
        _laminar_similarity_accuracy(api),
        check=check,
    )
    current &= write_or_check(
        ROOT / "docs/_generated/turbulent_sa_comparison.rst",
        _turbulent_sa_comparison(api),
        check=check,
    )
    current &= write_or_check(
        ROOT / "docs/_generated/marching_laminar_verification.rst",
        _marching_laminar_verification(api),
        check=check,
    )
    current &= write_or_check(TABLE_PATH, _rst(source, chart, invariants), check=check)
    current &= write_or_check(FRICTION_PATH, friction, check=check)
    current &= write_or_check(PROTRUSION_PATH, protrusion, check=check)
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not generate(check=args.check):
        raise SystemExit("viscous-flow verification artifacts are stale")


if __name__ == "__main__":
    main()
