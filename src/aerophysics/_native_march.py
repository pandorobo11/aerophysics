"""Optional native laminar Newton solve, edge evaluation, and downstream march."""

import importlib
import os
import warnings
from collections.abc import Callable
from typing import TYPE_CHECKING, Protocol, TypedDict, cast

import numpy as np

from aerophysics._array import FloatArray
from aerophysics.exceptions import ApplicabilityWarning, BoundaryLayerConvergenceError
from aerophysics.transport import DynamicViscosityModel, KeyesModel, SutherlandModel

if TYPE_CHECKING:
    from aerophysics._boundary_layer_march import (
        Discretization,
        MarchSolution,
        StationEdge,
    )
    from aerophysics.marching_boundary_layer import _EdgeInterpolant


class Kernel(Protocol):
    def __init__(
        self,
        eta: FloatArray,
        prandtl_number: float,
        viscosity_kind: int,
        viscosity_parameters: tuple[float, ...],
        jacobian: str = "analytic",
    ) -> None: ...

    def newton(
        self,
        xi: float,
        beta: float,
        heating: float,
        edge_temperature: float,
        wall_ratio: float | None,
        coefficient: float,
        lagged: FloatArray,
        guess: FloatArray,
        max_iterations: int,
        tolerance: float,
    ) -> tuple[FloatArray | None, bool]: ...


class EdgeKernel(Protocol):
    def states(self, distances: FloatArray) -> tuple[FloatArray, bool]: ...

    def polynomials(
        self, distances: FloatArray
    ) -> list[tuple[float, float, float | None]]: ...


class MarchOutput(TypedDict):
    status: str
    candidate: float
    positions: FloatArray
    states: FloatArray
    slopes: FloatArray
    edges: FloatArray
    accepted_positions: list[float]
    accepted_shear: list[float]
    out_of_range: bool


class Extension(Protocol):
    NativeEdge: Callable[..., EdgeKernel]
    march: Callable[..., MarchOutput]
    LaminarKernel: Callable[[FloatArray, float, int, tuple[float, ...], str], Kernel]


def load_extension() -> Extension | None:
    setting = os.environ.get("AEROPHYSICS_NATIVE", "auto")
    if setting not in ("", "auto", "0", "1"):
        raise ValueError("AEROPHYSICS_NATIVE must be auto, 0, or 1")
    if setting == "0":
        return None
    try:
        return cast(Extension, importlib.import_module("aerophysics_native"))
    except ImportError as error:
        if setting == "1":
            raise ImportError(
                "AEROPHYSICS_NATIVE=1 requires aerophysics_native"
            ) from error
        return None


def viscosity_spec(
    model: DynamicViscosityModel,
) -> tuple[int, tuple[float, ...]] | None:
    if type(model) is SutherlandModel:
        return 0, (
            model.reference_viscosity,
            model.reference_temperature,
            model.sutherland_temperature,
        )
    if type(model) is KeyesModel:
        return 1, (
            model.coefficient,
            model.additive_temperature,
            model.exponential_temperature,
            *model.temperature_range,
        )
    return None


class NativeMarch:
    def __init__(self, kernel: Kernel, model: DynamicViscosityModel) -> None:
        self.kernel = kernel
        self.model = model

    def newton(
        self,
        guess: FloatArray,
        edge: "StationEdge",
        coefficient: float,
        lagged: FloatArray,
        iterations: int,
        tolerance: float,
    ) -> FloatArray | None:
        state, out_of_range = self.kernel.newton(
            edge.xi,
            edge.beta,
            edge.heating,
            edge.temperature,
            edge.wall_ratio,
            coefficient,
            lagged,
            guess,
            iterations,
            tolerance,
        )
        self.warn_range(out_of_range)
        return state

    def warn_range(self, out_of_range: bool) -> None:
        if out_of_range and type(self.model) is KeyesModel:
            minimum, maximum = self.model.temperature_range
            warnings.warn(
                f"temperature is outside the fitted range {minimum:g}--"
                f"{maximum:g} K; the viscosity correlation was extrapolated",
                ApplicabilityWarning,
                stacklevel=3,
            )


def create_adapter(
    eta: FloatArray, prandtl_number: float, model: DynamicViscosityModel
) -> NativeMarch | None:
    extension = load_extension()
    spec = viscosity_spec(model)
    if extension is None or spec is None:
        return None
    jacobian = os.environ.get("AEROPHYSICS_NATIVE_JACOBIAN", "analytic") or "analytic"
    if jacobian not in ("analytic", "finite-difference"):
        raise ValueError(
            "AEROPHYSICS_NATIVE_JACOBIAN must be analytic or finite-difference"
        )
    kind, parameters = spec
    return NativeMarch(
        extension.LaminarKernel(eta, prandtl_number, kind, parameters, jacobian), model
    )


def create_edge(
    extension: Extension,
    interpolant: "_EdgeInterpolant",
    spec: tuple[int, tuple[float, ...]],
) -> EdgeKernel:
    from aerophysics.marching_boundary_layer import _GAUSS_NODES, _GAUSS_WEIGHTS

    pressure = interpolant.log_pressure
    slope = interpolant.log_pressure_slope
    wall = interpolant.wall_temperature
    gas = interpolant.gas
    return extension.NativeEdge(
        (pressure.x, pressure.c),
        (slope.x, slope.c),
        None if wall is None else (wall.x, wall.c),
        interpolant.panels,
        interpolant.xi_at_panels,
        _GAUSS_NODES,
        _GAUSS_WEIGHTS,
        (
            interpolant.total_temperature,
            interpolant.total_pressure,
            gas.heat_capacity_ratio,
            gas.cp,
            gas.specific_gas_constant,
        ),
        *spec,
    )


def native_march(
    interpolant: "_EdgeInterpolant",
    targets: FloatArray,
    discretization: "Discretization",
    model: DynamicViscosityModel,
) -> tuple["MarchSolution", FloatArray] | None:
    from aerophysics import _boundary_layer_march as reference

    setting = os.environ.get("AEROPHYSICS_NATIVE_MARCH", "1")
    if setting not in ("", "0", "1"):
        raise ValueError("AEROPHYSICS_NATIVE_MARCH must be 0 or 1")
    if setting == "0" or type(discretization) is not reference.Discretization:
        return None
    adapter = discretization._native
    spec = viscosity_spec(model)
    if adapter is None or spec is None:
        return None
    extension = load_extension()
    if extension is None:
        return None
    edge = create_edge(extension, interpolant, spec)
    raw = extension.march(
        adapter.kernel,
        edge,
        targets,
        discretization.eta,
        discretization.wall_weights,
        (
            reference._FIRST_STEP_FRACTION,
            reference._MINIMUM_STEP_FRACTION,
            reference._STEP_GROWTH,
            reference._MAX_STEP_RATIO,
            reference._MERGE_WITH_TARGET,
            reference._RELATIVE_STEP,
            reference._MAX_LOG_VELOCITY_CHANGE,
            reference._MAX_LOG_SHEAR_CHANGE,
            reference._MAX_HEAT_TRANSFER_CHANGE,
            reference._HEAT_TRANSFER_FLOOR,
            reference._SEPARATION_SHEAR_FRACTION,
            reference._EDGE_TOLERANCE,
        ),
        reference._EDGE_PROBE_NODES,
        discretization.newton_iterations,
        discretization.newton_tolerance,
    )
    adapter.warn_range(raw["out_of_range"])
    status = raw["status"]
    candidate = raw["candidate"]
    if status == "invalid_start":
        raise BoundaryLayerConvergenceError("the march must start at xi = 0")
    if status == "leading_edge_failed":
        raise BoundaryLayerConvergenceError(
            "the leading-edge similarity profile did not converge"
        )
    if status == "not_converged":
        raise BoundaryLayerConvergenceError(
            f"the march did not converge at distance {candidate:.6g}"
        )
    if status == "outer_boundary":
        raise BoundaryLayerConvergenceError(
            "the outer boundary is too close to the boundary layer at "
            f"distance {candidate:.6g}"
        )
    # Match _build_solution's field layout, including empty partial solutions.
    field = np.ascontiguousarray(raw["states"].transpose(1, 0, 2))
    solution = reference.MarchSolution(
        eta=discretization.eta,
        distance=raw["positions"],
        velocity_ratio=field[1],
        temperature_ratio=field[2],
        wall_velocity_slope=raw["slopes"][:, 0].copy(),
        wall_temperature_slope=raw["slopes"][:, 1].copy(),
    )
    if status == "separated":
        raise reference.SeparationError(
            reference._separation_position(
                raw["accepted_positions"], raw["accepted_shear"]
            ),
            solution,
        )
    return solution, raw["edges"]
