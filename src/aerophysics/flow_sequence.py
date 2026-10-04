"""Scalar sequences of local ideal-gas flow transformations.

All dimensional quantities use SI; angles are non-negative radians. Cone
steps hand off the cone *surface* state. No geometry or wave interaction is
solved. The gas and frozen composition are fixed throughout a sequence.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, log, sqrt
from typing import Literal

from scipy.optimize import brentq

from aerophysics._thermal_shocks import properties
from aerophysics.atmosphere import standard_atmosphere
from aerophysics.exceptions import (
    ModelRangeError,
    NoAttachedShockError,
)
from aerophysics.expansion import prandtl_meyer_expansion
from aerophysics.gas import AIR, PerfectGas
from aerophysics.isentropic import (
    MachBranch,
    area_ratio,
    isentropic_state,
    mach_from_area_ratio,
    mach_from_total_pressure_ratio,
)
from aerophysics.shocks import (
    NormalShockResult,
    ObliqueShockResult,
    ShockBranch,
    ShockGasModel,
    conical_shock,
    normal_shock,
    oblique_shock,
)


@dataclass(frozen=True, slots=True)
class FlowState:
    """Independent static inputs: pressure [Pa], temperature [K], Mach."""

    pressure: float
    temperature: float
    mach: float

    def __post_init__(self) -> None:
        for name, value in (
            ("pressure", self.pressure),
            ("temperature", self.temperature),
        ):
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not isfinite(self.mach) or self.mach < 0:
            raise ValueError("mach must be finite and non-negative")

    @classmethod
    def from_velocity(
        cls,
        pressure: float,
        temperature: float,
        velocity: float,
        gas: ShockGasModel = AIR,
    ) -> FlowState:
        """Initialize with speed magnitude [m/s]."""
        if not isfinite(velocity) or velocity < 0:
            raise ValueError("velocity must be finite and non-negative")
        return cls(pressure, temperature, velocity / _thermo(temperature, gas)[2])

    @classmethod
    def from_total(
        cls,
        total_pressure: float,
        total_temperature: float,
        mach: float,
        gas: ShockGasModel = AIR,
    ) -> FlowState:
        """Initialize with stagnation pressure [Pa] and temperature [K]."""
        state = isentropic_state(
            mach,
            gas,
            total_pressure=total_pressure,
            total_temperature=total_temperature,
            allow_extrapolation=False,
        )
        return cls(float(state.static_pressure), float(state.static_temperature), mach)

    @classmethod
    def from_atmosphere(cls, geometric_altitude: float, mach: float) -> FlowState:
        """Initialize p,T from USSA 1976 at geometric altitude [m]."""
        state = standard_atmosphere(geometric_altitude)
        return cls(float(state.pressure), float(state.temperature), mach)


@dataclass(frozen=True, slots=True)
class NormalShockStep:
    """Normal shock with no additional parameter."""

    name: str = "Normal shock"


@dataclass(frozen=True, slots=True)
class ObliqueShockStep:
    """Attached shock; deflection angle in radians."""

    angle: float
    branch: ShockBranch = ShockBranch.WEAK
    name: str = "Oblique shock"


@dataclass(frozen=True, slots=True)
class ConicalShockStep:
    """Weak attached cone solution; half-angle in radians."""

    angle: float
    name: str = "Cone surface"


@dataclass(frozen=True, slots=True)
class ExpansionStep:
    """Centered expansion; positive turn angle in radians."""

    angle: float
    name: str = "Expansion"


@dataclass(frozen=True, slots=True)
class IsentropicStep:
    """Target Mach, p2/p1, or A2/A1; area requires an explicit exit branch."""

    value: float
    basis: Literal["mach", "pressure_ratio", "area_ratio"] = "mach"
    branch: MachBranch | None = None
    name: str = "Isentropic"


type FlowStep = (
    NormalShockStep
    | ObliqueShockStep
    | ConicalShockStep
    | ExpansionStep
    | IsentropicStep
)


@dataclass(frozen=True, slots=True)
class FlowProperties:
    """Derived scalar SI properties; unavailable total conditions are None."""

    state: FlowState
    density: float
    velocity: float
    speed_of_sound: float
    dynamic_pressure: float
    mass_flux: float
    total_enthalpy: float
    total_temperature: float | None
    total_pressure: float | None
    warning: str = ""


@dataclass(frozen=True, slots=True)
class StepResult:
    """One exit state, or a failed/uncomputed stage with its explanation."""

    step: FlowStep
    status: str
    properties: FlowProperties | None = None
    pressure_ratio: float | None = None
    temperature_ratio: float | None = None
    total_pressure_recovery: float | None = None
    entropy_change: float | None = None
    shock_angle: float | None = None
    post_shock_mach: float | None = None
    message: str = ""


@dataclass(frozen=True, slots=True)
class FlowSequenceResult:
    """Initial properties and every requested stage, including failures."""

    initial: FlowProperties
    steps: tuple[StepResult, ...]

    @property
    def complete(self) -> bool:
        """Whether all requested stages succeeded (including an empty sequence)."""
        return all(step.status == "ok" for step in self.steps)


def _thermo(temperature: float, gas: ShockGasModel) -> tuple[float, float, float]:
    if not isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be finite and positive")
    if isinstance(gas, PerfectGas):
        return (
            gas.cp * temperature,
            gas.cp * log(temperature),
            float(gas.speed_of_sound(temperature)),
        )
    state = properties(temperature, gas)
    return state.enthalpy, state.entropy, sqrt(state.sound_speed_squared)


def flow_properties(state: FlowState, gas: ShockGasModel = AIR) -> FlowProperties:
    """Derive local and stagnation properties without thermal extrapolation."""
    h, s, sound = _thermo(state.temperature, gas)
    speed = state.mach * sound
    h0 = h + speed**2 / 2
    density = state.pressure / (gas.specific_gas_constant * state.temperature)
    t0: float | None
    p0: float | None
    warning = ""
    if isinstance(gas, PerfectGas):
        t0 = h0 / gas.cp
        p0 = state.pressure * (t0 / state.temperature) ** (
            gas.cp / gas.specific_gas_constant
        )
    elif h0 > _thermo(gas.temperature_range[1], gas)[0]:
        t0, p0 = None, None
        warning = (
            "Stagnation temperature exceeds the gas model range; "
            "total T and p unavailable"
        )
    else:
        t0 = float(
            brentq(
                lambda t: _thermo(t, gas)[0] - h0,
                state.temperature,
                gas.temperature_range[1],
                xtol=1e-10,
            )
        )
        p0 = state.pressure * exp((_thermo(t0, gas)[1] - s) / gas.specific_gas_constant)
    return FlowProperties(
        state,
        density,
        speed,
        sound,
        density * speed**2 / 2,
        density * speed,
        h0,
        t0,
        p0,
        warning,
    )


def _isentropic(
    state: FlowProperties, step: IsentropicStep, gas: ShockGasModel
) -> FlowState:
    if not isfinite(step.value) or step.value < 0:
        raise ValueError("isentropic target must be finite and non-negative")
    if state.total_temperature is None or state.total_pressure is None:
        raise ModelRangeError(
            "isentropic step requires stagnation conditions inside the model range"
        )
    t0, p0 = state.total_temperature, state.total_pressure
    if step.basis == "mach":
        mach = step.value
    elif step.basis == "pressure_ratio":
        if step.value <= 0:
            raise ValueError("pressure ratio must be positive")
        mach = float(
            mach_from_total_pressure_ratio(
                p0 / (state.state.pressure * step.value),
                gas,
                total_temperature=t0,
                allow_extrapolation=False,
            )
        )
    elif step.basis == "area_ratio":
        if step.value <= 0 or step.branch is None:
            raise ValueError("area ratio must be positive and requires an exit branch")
        target = step.value * float(
            area_ratio(
                state.state.mach, gas, total_temperature=t0, allow_extrapolation=False
            )
        )
        mach = float(
            mach_from_area_ratio(
                target,
                MachBranch(step.branch),
                gas,
                total_temperature=t0,
                allow_extrapolation=False,
            )
        )
    else:
        raise ValueError("unknown isentropic basis")
    return FlowState.from_total(p0, t0, mach, gas)


def _advance(
    state: FlowProperties, step: FlowStep, gas: ShockGasModel
) -> tuple[FlowState, float, float | None, float | None]:
    upstream = state.state
    if isinstance(step, IsentropicStep):
        return _isentropic(state, step, gas), 1.0, None, None
    if not isinstance(
        step, (NormalShockStep, ObliqueShockStep, ConicalShockStep, ExpansionStep)
    ):
        raise TypeError("unsupported flow step")
    if upstream.mach < 1:
        raise ValueError("shock and expansion steps require upstream Mach >= 1")
    if not isinstance(step, NormalShockStep) and (
        not isfinite(step.angle) or step.angle < 0
    ):
        raise ValueError("angle must be finite and non-negative radians")
    if isinstance(step, ConicalShockStep):
        cone = conical_shock(
            upstream.mach, step.angle, gas, upstream_temperature=upstream.temperature
        )
        return (
            FlowState(
                upstream.pressure * float(cone.surface_pressure_ratio),
                upstream.temperature * float(cone.surface_temperature_ratio),
                float(cone.surface_mach),
            ),
            float(cone.total_pressure_ratio),
            float(cone.shock_angle),
            float(cone.post_shock_mach),
        )
    angle = None
    shock: NormalShockResult | ObliqueShockResult
    if isinstance(step, NormalShockStep):
        shock = normal_shock(
            upstream.mach, gas, upstream_temperature=upstream.temperature
        )
    elif isinstance(step, ObliqueShockStep):
        oblique = oblique_shock(
            upstream.mach,
            step.angle,
            ShockBranch(step.branch),
            gas,
            upstream_temperature=upstream.temperature,
        )
        angle = float(oblique.shock_angle)
        shock = oblique
    else:
        expansion = prandtl_meyer_expansion(
            upstream.mach, step.angle, gas, upstream_temperature=upstream.temperature
        )
        return (
            FlowState(
                upstream.pressure * float(expansion.static_pressure_ratio),
                upstream.temperature * float(expansion.static_temperature_ratio),
                float(expansion.downstream_mach),
            ),
            1.0,
            None,
            None,
        )
    return (
        FlowState(
            upstream.pressure * float(shock.static_pressure_ratio),
            upstream.temperature * float(shock.static_temperature_ratio),
            float(shock.downstream_mach),
        ),
        float(shock.total_pressure_ratio),
        angle,
        None,
    )


def solve_flow_sequence(
    initial: FlowState,
    steps: list[FlowStep] | tuple[FlowStep, ...],
    gas: ShockGasModel = AIR,
) -> FlowSequenceResult:
    """Apply stages in order; preserve valid prefixes and classify stage failures.

    Invalid initial states raise. A failed stage stops propagation; all later
    stages have status ``not_computed``. Entropy change [J/(kg K)] and total
    pressure recovery are relative to the initial state.
    """
    first = current = flow_properties(initial, gas)
    results: list[StepResult] = []
    recovery = 1.0
    entropy = 0.0
    failed = False
    for step in steps:
        if failed:
            results.append(StepResult(step, "not_computed"))
            continue
        try:
            output, loss, angle, post_mach = _advance(current, step, gas)
            derived = flow_properties(output, gas)
            recovery *= loss
            entropy -= gas.specific_gas_constant * log(loss)
            results.append(
                StepResult(
                    step,
                    "ok",
                    derived,
                    output.pressure / current.state.pressure,
                    output.temperature / current.state.temperature,
                    recovery,
                    entropy,
                    angle,
                    post_mach,
                )
            )
            current = derived
        except (ValueError, RuntimeError) as error:
            status = (
                "detached"
                if isinstance(error, NoAttachedShockError)
                else "out_of_range"
                if isinstance(error, ModelRangeError)
                else "convergence_error"
                if isinstance(error, RuntimeError)
                else "invalid_input"
            )
            results.append(StepResult(step, status, message=str(error)))
            failed = True
    return FlowSequenceResult(first, tuple(results))
