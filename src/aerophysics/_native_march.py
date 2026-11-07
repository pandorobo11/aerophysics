"""Optional native implementation of the laminar station Newton solve."""

import importlib
import os
import warnings
from collections.abc import Callable
from typing import TYPE_CHECKING, Protocol, cast

from aerophysics._array import FloatArray
from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.transport import DynamicViscosityModel, KeyesModel, SutherlandModel

if TYPE_CHECKING:
    from aerophysics._boundary_layer_march import StationEdge


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


class Extension(Protocol):
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
        if out_of_range and type(self.model) is KeyesModel:
            minimum, maximum = self.model.temperature_range
            warnings.warn(
                f"temperature is outside the fitted range {minimum:g}--"
                f"{maximum:g} K; the viscosity correlation was extrapolated",
                ApplicabilityWarning,
                stacklevel=2,
            )
        return state


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
