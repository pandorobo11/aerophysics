"""Generate boundary-layer and protrusion verification artifacts."""

from __future__ import annotations

import argparse
import csv
import sys
import warnings
from itertools import pairwise
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
    from aerophysics._boundary_layer_march import (
        Discretization,
        StationEdge,
        eta_grid,
        march,
    )
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
            Discretization(eta, 1.0, _LinearViscosity()),
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
            Discretization(eta, 0.72, wall_edge.viscosity_model),
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


def _marching_sa_verification() -> str:
    """Frozen ZPG regression and an unfitted external same-model comparison."""
    from aerophysics import (
        AIR,
        BoundaryLayerEdge,
        BoundaryLayerRegime,
        marching_boundary_layer,
    )
    from aerophysics.transport import AIR_VISCOSITY

    lines = [
        ".. list-table:: SA-noft2 ZPG differences from the frozen reference [%]",
        "   :header-rows: 1",
        "",
        "   * - Case",
        "     - Re_s",
        "     - Cf",
        "     - delta99",
        "     - theta",
        "     - Heat coefficient",
    ]
    with (REFERENCE / "turbulent_sa_boundary_layer.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    for label in dict.fromkeys(row["case"] for row in rows):
        group = [row for row in rows if row["case"] == label]
        first = group[0]
        mach, temperature = float(first["mach"]), float(first["edge_temperature_K"])
        velocity = mach * np.sqrt(1.4 * AIR.specific_gas_constant * temperature)
        density = 1e7 * AIR_VISCOSITY.dynamic_viscosity(temperature) / velocity
        pressure = density * AIR.specific_gas_constant * temperature
        distances = np.array([float(row["reynolds_number"]) / 1e7 for row in group])
        ratio = 1.0 + 0.2 * mach**2
        edge = BoundaryLayerEdge.from_pressure(
            np.r_[0.0, distances],
            np.full(len(distances) + 1, pressure),
            total_pressure=pressure * ratio**3.5,
            total_temperature=temperature * ratio,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = marching_boundary_layer(
                edge,
                regime=BoundaryLayerRegime.TURBULENT,
                wall_temperature=float(first["wall_temperature_K"]),
            )
        quantities = (
            result.local_skin_friction_coefficient,
            result.boundary_layer_thickness / distances,
            result.momentum_thickness / distances,
            result.wall_heat_flux / (density * velocity * AIR.cp * temperature),
        )
        keys = ("cf", "delta99_over_x", "momentum_over_x", "wall_heat_flux_coefficient")
        for index, row in enumerate(group):
            differences = []
            for quantity, key in zip(quantities, keys, strict=True):
                reference = float(row[key])
                scale = (
                    float(row["cf"])
                    if key == keys[-1] and "AD" in label
                    else abs(reference)
                )
                difference = 100.0 * (quantity[index] - reference) / scale
                if abs(difference) >= 0.3:
                    raise AssertionError(
                        f"SA ZPG {label} Re={row['reynolds_number']} {key}: "
                        f"{difference:.6g}% exceeds 0.3%"
                    )
                differences.append(difference)
            lines += [f"   * - {label}", f"     - {float(row['reynolds_number']):.0e}"]
            lines += [f"     - {value:.4f}" for value in differences]
    lines += [
        "",
        "The recovery-temperature rows (AD) use an absolute heat-coefficient",
        "difference divided by the reference Cf. All other differences are relative",
        "to the reference quantity. The limit is 0.3% per entry; reference grid",
        "uncertainty is about 0.1%. The reference displacement integrates only to",
        "delta99 and is excluded here. Unit Reynolds number is 1e7 per metre,",
        "gamma is 1.4, Pr is 0.72, Prt is 0.9, and Sutherland S is 110.4 K.",
        "",
        ".. list-table:: NASA TMR bump: local wall shear / freestream dynamic pressure",
        "   :header-rows: 1",
        "",
        "   * - Cp source / Cf reference",
        "     - Region (x)",
        "     - Median difference [%]",
        "     - Maximum absolute difference [%]",
        "     - Maximum code spread [%]",
    ]
    for row in _tmr_bump_comparison():
        if row[3] >= 3.0:
            raise AssertionError(
                f"SA TMR {row[0]} {row[1]}: {row[3]:.6g}% exceeds 3.0%"
            )
        lines += [f"   * - {row[0]}", f"     - {row[1]}"]
        lines += [f"     - {value:.3f}" for value in row[2:]]
    lines += [
        "",
        "Each march uses that code's wall Cp as edge Cp. The comparison excludes",
        "x < 0.05 and x >= 1.45; the march ends at the last source point below 1.48.",
        "Cp is replaced by its linearly interpolated value at x = 0.02 for x < 0.02,",
        "including the added leading edge x = 0. No other smoothing or fitting",
        "is used.",
        "The wall distance is the trapezoidal integral of sqrt(1 + (dy/dx)^2) on",
        "the source x nodes. Regions are upstream [0.05, 0.3), windward [0.3, 0.675),",
        "crest [0.675, 0.825), leeward [0.825, 1.2), recovery [1.2, 1.45).",
        "",
        "M_inf = 0.2 and Re_inf = 3e6 per metre. Dimensional edge nu_tilde is",
        "3 nu_inf, converted to the initial local edge ratio. T_inf = 300 K,",
        "an adiabatic wall,",
        "and freestream dynamic-pressure normalization for both Cp and Cf are",
        "assumptions: the source does not state these three details. Repeating at",
        "250 K and 350 K changes each regional maximum by less than 0.007 percentage",
        "points; temperature has little effect at this Mach number. The small",
        "upstream-flat discrepancy supports the normalization assumption. Wall Cp",
        "need not equal boundary-layer edge Cp on the curved bump. The regression",
        "limit is 3.0%, from the measured maximum of 2.50% plus 0.50 percentage",
        "points of margin; it is not an accuracy claim. The same-model CFD spread",
        "is smaller than our differences over much of the bump. Accuracy against",
        "experiment is not established by this comparison.",
        "",
    ]
    return "\n".join(lines)


def _tmr_bump_comparison(
    temperature: float = 300.0,
) -> list[tuple[str, str, float, float, float]]:
    from aerophysics import (
        AIR,
        BoundaryLayerEdge,
        BoundaryLayerRegime,
        marching_boundary_layer,
    )
    from aerophysics.transport import AIR_VISCOSITY

    def zones(name: str) -> dict[str, np.ndarray]:
        records: dict[str, list[list[float]]] = {}
        label = ""
        for line in (REFERENCE / name).read_text().splitlines():
            if line.lower().startswith("zone"):
                label = line.split('"')[1]
                records[label] = []
            elif label:
                records[label].append([float(value) for value in line.split()])
        return {key: np.array(value) for key, value in records.items()}

    cp, cf = zones("tmr_bump_sa_cp.dat"), zones("tmr_bump_sa_cf.dat")
    rows = []
    for code in cp:
        x = np.unique(np.r_[0.0, cp[code][:, 0]])
        x = x[(x >= 0.0) & (x <= 1.48)]
        coefficient = np.interp(np.maximum(x, 0.02), cp[code][:, 0], cp[code][:, 1])
        angle = np.pi * x / 0.9 - np.pi / 3.0
        slope = np.where(
            (x >= 0.3) & (x <= 1.2),
            0.2 * np.pi / 0.9 * np.sin(angle) ** 3 * np.cos(angle),
            0.0,
        )
        distance = cumulative_trapezoid(np.sqrt(1.0 + slope**2), x, initial=0.0)
        mach = 0.2
        velocity = mach * np.sqrt(1.4 * AIR.specific_gas_constant * temperature)
        density = 3e6 * AIR_VISCOSITY.dynamic_viscosity(temperature) / velocity
        pressure = density * AIR.specific_gas_constant * temperature
        dynamic = 0.5 * density * velocity**2
        edge = BoundaryLayerEdge.from_pressure(
            distance,
            pressure + dynamic * coefficient,
            total_pressure=pressure * (1.0 + 0.2 * mach**2) ** 3.5,
            total_temperature=temperature * (1.0 + 0.2 * mach**2),
        )
        initial_nu = (
            AIR_VISCOSITY.dynamic_viscosity(float(edge.temperature[0]))
            * AIR.specific_gas_constant
            * edge.temperature[0]
            / edge.pressure[0]
        )
        farfield_nu = AIR_VISCOSITY.dynamic_viscosity(temperature) / density
        result = marching_boundary_layer(
            edge,
            regime=BoundaryLayerRegime.TURBULENT,
            edge_sa_ratio=3.0 * farfield_nu / initial_nu,
        )
        computed = result.wall_shear_stress / dynamic
        reference = np.interp(x[1:], cf[code][:, 0], cf[code][:, 1])
        other_code = "FUN3D" if code == "CFL3D" else "CFL3D"
        other = np.interp(x[1:], cf[other_code][:, 0], cf[other_code][:, 1])
        difference, spread = (
            100.0 * (computed / reference - 1.0),
            100.0 * (other / reference - 1.0),
        )
        for region, left, right in (
            ("Upstream", 0.05, 0.3),
            ("Windward", 0.3, 0.675),
            ("Crest", 0.675, 0.825),
            ("Leeward", 0.825, 1.2),
            ("Recovery", 1.2, 1.45),
        ):
            mask = (x[1:] >= left) & (x[1:] < right)
            rows.append(
                (
                    code,
                    f"{region} [{left:g}, {right:g})",
                    float(np.median(difference[mask])),
                    float(np.max(abs(difference[mask]))),
                    float(np.max(abs(spread[mask]))),
                )
            )
    return rows


# R&M 3633 (Winter, Rotta and Smith, 1968), Fig. 1. Coefficients in
# descending powers, eta=R/l and x=X/l, l=60 in. The last segment is CUBIC.
_RM_JUNCTIONS = np.array([0.0, 0.14158, 0.45725, 0.64473, 0.76317, 1.0])
_RM_POLYNOMIALS = (
    (0.36397, 0.0),
    (6.34730, -6.49221, 1.22668, 0.33498, -0.00461),
    (-0.236382, 0.607620, -0.595521, 0.229694),
    (-53.58258, 150.87806, -158.03866, 72.96859, -12.4939),
    (0.236382, -0.390789, 0.290242, -0.064001),
)


def _rm_geometry(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Radius and accurately integrated wall distance, in metres."""
    from scipy.integrate import quad

    gaps = [
        abs(np.polyval(a, t) - np.polyval(b, t))
        for a, b, t in zip(
            _RM_POLYNOMIALS[:-1], _RM_POLYNOMIALS[1:], _RM_JUNCTIONS[1:-1], strict=True
        )
    ]
    if max(gaps) >= 1e-5:
        raise AssertionError(f"R&M 3633 geometry transcription: junction gaps {gaps}")
    polynomials = [np.poly1d(c) for c in _RM_POLYNOMIALS]
    slopes = [p.deriv() for p in polynomials]

    def integral(a: float, b: float, i: int) -> float:
        return float(
            quad(
                lambda t: np.sqrt(1 + slopes[i](t) ** 2),
                a,
                b,
                epsabs=1e-12,
                epsrel=1e-12,
            )[0]
        )

    lengths = np.r_[
        0.0,
        np.cumsum(
            [integral(a, b, i) for i, (a, b) in enumerate(pairwise(_RM_JUNCTIONS))]
        ),
    ]
    indices = np.minimum(np.searchsorted(_RM_JUNCTIONS, x, side="right") - 1, 4)
    radius = (
        np.array([polynomials[i](t) for i, t in zip(indices, x, strict=True)]) * 1.524
    )
    distance = (
        np.array(
            [
                lengths[i] + integral(float(_RM_JUNCTIONS[i]), float(t), int(i))
                for i, t in zip(indices, x, strict=True)
            ]
        )
        * 1.524
    )
    return radius, distance


def _rm3633_comparison(
    total_temperature: float = 310.0, *, cone_pressure: bool = False
) -> dict:
    """Unfitted wind-tunnel comparison; retain failures and applicability warnings."""
    from scipy.interpolate import PchipInterpolator

    from aerophysics import (
        AIR,
        BoundaryLayerEdge,
        BoundaryLayerRegime,
        conical_shock,
        marching_boundary_layer,
    )
    from aerophysics.exceptions import (
        BoundaryLayerConvergenceError,
        BoundaryLayerSeparationError,
    )
    from aerophysics.transport import AIR_VISCOSITY

    measured = np.genfromtxt(REFERENCE / "rm3633_m17.csv", delimiter=",")
    if measured.shape != (28, 4):
        raise AssertionError("R&M 3633 transcription must have 28 rows and 4 columns")
    total_ratios = (
        (1 + 0.7 * 1.7**2 * measured[:, 1])
        * (1 + 0.2 * measured[:, 2] ** 2) ** 3.5
        / (1 + 0.2 * 1.7**2) ** 3.5
    )
    if not np.all((total_ratios > 0.9925) & (total_ratios < 0.9963)):
        raise AssertionError("R&M 3633 Cp/M_delta transcription consistency failed")
    cone = conical_shock(1.7, np.radians(20.0))
    cone_cp = (float(cone.surface_pressure_ratio) - 1.0) / (0.7 * 1.7**2)
    temperature = total_temperature / (1 + 0.2 * 1.7**2)
    velocity = 1.7 * np.sqrt(1.4 * AIR.specific_gas_constant * temperature)
    density = 10.08e6 / 1.524 * AIR_VISCOSITY.dynamic_viscosity(temperature) / velocity
    pressure = density * AIR.specific_gas_constant * temperature
    # Preserve measured values, including the unsmoothed small variations.
    source_x = np.r_[0.0, measured[:, 0] / 60.0]
    cp = np.r_[measured[0, 1], measured[:, 1]]
    if cone_pressure:
        cp[source_x <= 0.1] = cone_cp
    _, source_s = _rm_geometry(source_x)
    log_p = PchipInterpolator(source_s, np.log(pressure * (1 + 0.7 * 1.7**2 * cp)))
    # About 0.25 in spacing, including every measured station and body junction.
    x = np.unique(np.r_[np.arange(237) / 240.0, source_x, _RM_JUNCTIONS[1:-1]])
    radius, distance = _rm_geometry(x)
    edge = BoundaryLayerEdge.from_pressure(
        distance,
        np.exp(log_p(distance)),
        total_pressure=pressure
        * (1 + 0.2 * 1.7**2) ** 3.5
        * float(cone.total_pressure_ratio),
        total_temperature=total_temperature,
    )
    failure = None
    result = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = marching_boundary_layer(
                edge,
                regime=BoundaryLayerRegime.TURBULENT,
                wall_radius=radius,
                edge_sa_ratio=3.0,
            )
        except BoundaryLayerSeparationError as error:
            failure = str(error)
            result = error.result
        except BoundaryLayerConvergenceError as error:
            failure = str(error)
    record = {
        "temperature": total_temperature,
        "cone_pressure": cone_pressure,
        "cone_cp": cone_cp,
        "cone_total_pressure_ratio": float(cone.total_pressure_ratio),
        "freestream_pressure_Pa": float(pressure),
        "failure": failure,
        "warnings": [str(w.message) for w in caught],
        "warning_kinds": sorted(
            {
                f"{w.category.__name__}: transverse curvature is neglected"
                if "transverse curvature is neglected" in str(w.message)
                else w.category.__name__
                for w in caught
            }
        ),
        "rows": [],
    }
    if result is not None and result.distance.size:
        available = (
            (measured[:, 0] >= 12)
            & np.isfinite(measured[:, 3])
            & (source_s[1:] <= result.distance[-1])
        )
        stations = measured[available]
        computed = np.interp(
            source_s[1:][available],
            result.distance,
            result.local_skin_friction_coefficient,
        )
        ratios = np.interp(result.distance, distance, radius)
        curvature = result.boundary_layer_thickness / ratios
        record["max_delta_over_radius"] = float(np.max(curvature))
        record["max_delta_over_radius_s_m"] = float(
            result.distance[np.argmax(curvature)]
        )
        record["curvature_at_stations"] = np.interp(
            source_s[1:][available], result.distance, curvature
        ).tolist()
        record["rows"] = [
            (
                float(row[0]),
                float(row[3] / 1000),
                float(value),
                float(value - row[3] / 1000),
            )
            for row, value in zip(stations, computed, strict=True)
        ]
        if computed.size:
            record["max_absolute_difference"] = float(
                np.max(abs(computed - stations[:, 3] / 1000))
            )
    return record


def _axisymmetric_display(value: float, digits: int, *, signed: bool = False) -> str:
    """Reject values within 1e-3 display units of a rounding boundary."""
    scaled = value * 10**digits
    distance = abs(scaled - np.floor(scaled) - 0.5)
    if not np.isfinite(scaled) or distance < 1e-3:
        raise AssertionError(
            f"axisymmetric display rounding boundary: value={value!r}, "
            f"digits={digits}, distance={distance!r} display units"
        )
    return format(value, f"{'+' if signed else ''}.{digits}f")


def _marching_axisymmetric_verification() -> str:
    """Generated V19 table; geometry and comparison regressions run on regeneration."""
    lines = [
        ".. list-table:: R&M 3633, M_inf=1.7, Re_l=10.08e6: local Cf",
        "   :header-rows: 1",
        "",
        "   * - X [in]",
        "     - Measured 10^3 cf",
        "     - SA-noft2 10^3 cf",
        "     - Difference (10^3 cf)",
        "     - abs difference / 0.0002",
        "     - delta99/r_w",
    ]
    base = _rm3633_comparison()
    for row, curvature in zip(
        base["rows"], base.get("curvature_at_stations", []), strict=True
    ):
        x, measured, computed, difference = row
        lines += [
            f"   * - {_axisymmetric_display(x, 1).removesuffix('.0')}",
            f"     - {_axisymmetric_display(1000 * measured, 2)}",
            f"     - {_axisymmetric_display(1000 * computed, 2)}",
            f"     - {_axisymmetric_display(1000 * difference, 2, signed=True)}",
            f"     - {_axisymmetric_display(abs(difference) / 0.0002, 1)}",
            f"     - {_axisymmetric_display(curvature, 3)}",
        ]
    lines += [
        "",
        "Measured Cf accuracy is +/-0.0002; Cp accuracy is +/-0.005 (section 3.6).",
        "The X=15 in Cf is missing. X<12 in is excluded because of transition.",
        "This is a wind-tunnel comparison, not a calibrated accuracy claim.",
        "",
        ".. list-table:: R&M 3633 sensitivity to assumptions",
        "   :header-rows: 1",
        "",
        "   * - Assumption",
        "     - Maximum absolute difference (10^3 cf)",
        "     - Maximum delta99/r_w",
        "     - Outcome",
    ]
    for label, record in [
        ("T0=310 K", base),
        ("T0=290 K", _rm3633_comparison(290.0)),
        ("T0=330 K", _rm3633_comparison(330.0)),
        ("T0=310 K, cone Cp replaced", _rm3633_comparison(cone_pressure=True)),
    ]:
        # This is a regression bound, not the measurement-accuracy tolerance.
        if record["failure"] is not None or record["max_absolute_difference"] > 0.0006:
            raise AssertionError(f"R&M 3633 regression: {record}")
        lines += [
            f"   * - {label}",
            "     - "
            + _axisymmetric_display(1000 * record["max_absolute_difference"], 2),
            f"     - {_axisymmetric_display(record['max_delta_over_radius'], 3)}",
            f"     - {record['failure'] or 'Completed'}",
        ]
    lines += ["", "Warnings from the baseline:", ""]
    lines += [
        f"* {kind}; maximum delta99/r_w = "
        f"{_axisymmetric_display(base['max_delta_over_radius'], 3)} "
        f"at s = {_axisymmetric_display(base['max_delta_over_radius_s_m'], 3)} m."
        for kind in base["warning_kinds"]
    ] or ["* None."]
    lines += [
        "",
        "The regression bound is 0.0006: the measured maximum over these cases is",
        "rounded up to 0.0005, plus 0.0001 margin. This is not a physical accuracy",
        "bound. Input construction and limitations are in the boundary-layer guide.",
        "",
    ]
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
    current &= write_or_check(
        ROOT / "docs/_generated/marching_sa_verification.rst",
        _marching_sa_verification(),
        check=check,
    )
    current &= write_or_check(
        ROOT / "docs/_generated/marching_axisymmetric_verification.rst",
        _marching_axisymmetric_verification(),
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
