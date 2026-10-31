"""Laminar sphere and laminar/turbulent swept-cylinder heating in SI units.

References
----------
Tauber, M. E., *A Review of High-Speed, Convective, Heat-Transfer
Computation Methods*, NASA TP-2914, 1989.
https://ntrs.nasa.gov/citations/19890017745
Fay, J. A. and Riddell, F. R., *Theory of Stagnation Point Heat Transfer in
Dissociated Air*, Journal of the Aeronautical Sciences, 25(2), 73--85, 1958.
West, T. K., IV and Brandis, A. M., *Updated Stagnation Point Aeroheating
Correlations for Mars Entry*, NASA NTRS document 20200002354.
https://ntrs.nasa.gov/citations/20200002354
Holden, M. S., Moselle, J. R., and Lee, J., *Studies of Aerothermal Loads
Generated in Regions of Shock/Shock Interaction in Hypersonic Flow*, NASA
CR-181893, 1991, equations (5.7)--(5.12).
https://ntrs.nasa.gov/citations/19920002101
"""

import warnings
from dataclasses import dataclass
from enum import StrEnum
from typing import Final

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.gas import AIR, PerfectGas
from aerophysics.shocks import normal_shock
from aerophysics.transport import AIR_VISCOSITY, DynamicViscosityModel


class StagnationHeatingModel(StrEnum):
    """Sphere stagnation-point convective heating correlation."""

    FAY_RIDDELL = "fay_riddell"
    TAUBER = "tauber"


class AttachmentLineHeatingModel(StrEnum):
    """Turbulent swept-cylinder attachment-line heating correlation."""

    POLL = "poll"  # NASA CR-181893 equations (5.9)--(5.10)


POLL_TURBULENT_COEFFICIENT: Final = 0.0689
POLL_TURBULENT_EXPONENT: Final = 0.42
DEFAULT_ATTACHMENT_LINE_RECOVERY_FACTOR: Final = 0.89
BENARD_MINIMUM_EDGE_MACH: Final = 0.35
BENARD_MAXIMUM_EDGE_MACH: Final = 6.89
BENARD_MINIMUM_WALL_TOTAL_RATIO: Final = 0.14
BENARD_MAXIMUM_WALL_TOTAL_RATIO: Final = 0.96
FAY_RIDDELL_SPHERE_COEFFICIENT: Final = 0.763
TAUBER_SPHERE_COEFFICIENT: Final = 1.83e-4
TAUBER_CYLINDER_COEFFICIENT: Final = 1.29e-4
TAUBER_SWEEP_FACTOR: Final = 0.18
PERFECT_GAS_STAGNATION_TEMPERATURE_LIMIT: Final = 2000.0
TAUBER_MAXIMUM_VELOCITY: Final = 14_000.0
TAUBER_MAXIMUM_WALL_TEMPERATURE: Final = 4000.0
# h_inf / (V_inf^2/2), TP-2914 high-speed assumption
TAUBER_MAXIMUM_ENTHALPY_RATIO: Final = 0.1


@dataclass(frozen=True, slots=True)
class StagnationHeatingResult:
    """Sphere heating and freestream/stagnation state in SI units.

    Heat flux is in W/m², positive into the wall; velocity gradient is in 1/s.
    """

    model: StagnationHeatingModel
    heat_flux: FloatResult
    freestream_velocity: FloatResult
    freestream_density: FloatResult
    stagnation_temperature: FloatResult
    stagnation_pressure: FloatResult
    stagnation_enthalpy: FloatResult
    wall_temperature: FloatResult
    velocity_gradient: FloatResult | None


@dataclass(frozen=True, slots=True)
class SweptCylinderHeatingResult:
    """Laminar cylinder heating and flow state in SI units.

    Heat flux is in W/m², positive into the wall; sweep angle is in radians.
    """

    heat_flux: FloatResult
    sweep_angle: FloatResult
    freestream_velocity: FloatResult
    freestream_density: FloatResult
    recovery_enthalpy: FloatResult
    wall_temperature: FloatResult


@dataclass(frozen=True, slots=True)
class AttachmentLineHeatingResult:
    """Turbulent attachment-line heating and edge state in SI units.

    Heat flux is in W/m², positive into the wall; sweep angle is in radians.
    ``stanton_number`` is q_w / [rho* W_e c_p (T_aw - T_w)] with rho* the
    density at the edge pressure and reference temperature.
    """

    model: AttachmentLineHeatingModel
    heat_flux: FloatResult
    sweep_angle: FloatResult
    freestream_velocity: FloatResult
    freestream_density: FloatResult
    edge_temperature: FloatResult
    edge_pressure: FloatResult
    edge_density: FloatResult
    edge_velocity: FloatResult
    edge_mach: FloatResult
    velocity_gradient: FloatResult
    recovery_temperature: FloatResult
    reference_temperature: FloatResult
    reynolds_number: FloatResult
    stanton_number: FloatResult
    wall_temperature: FloatResult


def _inputs(*values: tuple[str, ArrayLike]) -> tuple[list[FloatArray], bool]:
    arrays = []
    scalar = True
    for name, value in values:
        array, is_scalar = as_float_array(value, name=name)
        if name == "mach":
            if np.any(array <= 1):
                raise ValueError("mach must be greater than one")
        elif name != "sweep_angle" and np.any(array <= 0):
            raise ValueError(f"{name} must be greater than zero")
        arrays.append(array)
        scalar = scalar and is_scalar
    try:
        return list(np.broadcast_arrays(*arrays)), scalar
    except ValueError as error:
        raise ValueError("heating inputs must be broadcastable") from error


def _freestream(
    mach: FloatArray, temperature: FloatArray, pressure: FloatArray, gas: PerfectGas
) -> tuple[FloatArray, FloatArray]:
    density = pressure / (gas.specific_gas_constant * temperature)
    velocity = mach * np.sqrt(
        gas.heat_capacity_ratio * gas.specific_gas_constant * temperature
    )
    return density, velocity


def _tauber_warnings(
    velocity: FloatArray, wall: FloatArray, freestream_enthalpy: FloatArray
) -> None:
    if np.any(velocity > TAUBER_MAXIMUM_VELOCITY):
        warnings.warn(
            "Tauber correlation exceeds the 14 km/s data range",
            ApplicabilityWarning,
            stacklevel=3,
        )
    if np.any(wall > TAUBER_MAXIMUM_WALL_TEMPERATURE):
        warnings.warn(
            "Tauber correlation exceeds the 4000 K wall-temperature assumption",
            ApplicabilityWarning,
            stacklevel=3,
        )
    if np.any(
        freestream_enthalpy / (0.5 * velocity**2) > TAUBER_MAXIMUM_ENTHALPY_RATIO
    ):
        warnings.warn(
            "Tauber correlation assumes h_inf << V^2/2; freestream enthalpy "
            "exceeds 0.1 of the kinetic energy",
            ApplicabilityWarning,
            stacklevel=3,
        )


def sphere_stagnation_heat_flux(
    mach: ArrayLike,
    temperature: ArrayLike,
    pressure: ArrayLike,
    nose_radius: ArrayLike,
    *,
    wall_temperature: ArrayLike,
    model: StagnationHeatingModel = StagnationHeatingModel.FAY_RIDDELL,
    prandtl_number: float = 0.72,
    gas: PerfectGas = AIR,
    viscosity_model: DynamicViscosityModel = AIR_VISCOSITY,
) -> StagnationHeatingResult:
    """Estimate laminar sphere stagnation heat flux, positive into the wall.

    Inputs are SI (temperatures in K, pressure in Pa, radius in m), with
    supersonic freestream Mach number. Fay--Riddell assumes Lewis number one,
    no dissociation and a Newtonian velocity gradient. Tauber's equation 40
    is an air correlation; ``gas`` only converts the freestream state and
    ``prandtl_number`` and ``viscosity_model`` are ignored for that model.
    See :doc:`../models/aerodynamic_heating` for applicability limits.
    """
    if not isinstance(model, StagnationHeatingModel):
        raise ValueError("model must be a StagnationHeatingModel value")
    if model is StagnationHeatingModel.FAY_RIDDELL and (
        not np.isfinite(prandtl_number) or prandtl_number <= 0
    ):
        raise ValueError("prandtl_number must be finite and greater than zero")
    (m, t, p, radius, wall), scalar = _inputs(
        ("mach", mach),
        ("temperature", temperature),
        ("pressure", pressure),
        ("nose_radius", nose_radius),
        ("wall_temperature", wall_temperature),
    )
    density, velocity = _freestream(m, t, p, gas)
    gamma = gas.heat_capacity_ratio
    factor = 1 + (gamma - 1) / 2 * m**2
    t0 = t * factor
    ps = (
        p
        * factor ** (gamma / (gamma - 1))
        * np.asarray(normal_shock(m, gas).total_pressure_ratio)
    )
    enthalpy = gas.cp * t0
    gradient = None
    if model is StagnationHeatingModel.FAY_RIDDELL:
        if np.any(t0 > PERFECT_GAS_STAGNATION_TEMPERATURE_LIMIT):
            warnings.warn(
                "Fay-Riddell edge properties use a perfect gas above 2000 K "
                "stagnation temperature; dissociation is not modelled",
                ApplicabilityWarning,
                stacklevel=2,
            )
        rho_e = ps / (gas.specific_gas_constant * t0)
        rho_w = ps / (gas.specific_gas_constant * wall)
        mu_e = viscosity_model.dynamic_viscosity(t0)
        mu_w = viscosity_model.dynamic_viscosity(wall)
        grad = np.sqrt(2 * (ps - p) / rho_e) / radius
        flux = (
            FAY_RIDDELL_SPHERE_COEFFICIENT
            * prandtl_number**-0.6
            * (rho_e * mu_e) ** 0.4
            * (rho_w * mu_w) ** 0.1
            * np.sqrt(grad)
            * gas.cp
            * (t0 - wall)
        )
        gradient = return_float(grad, scalar=scalar)
    else:
        _tauber_warnings(velocity, wall, gas.cp * t)
        flux = (
            TAUBER_SPHERE_COEFFICIENT
            * np.sqrt(density / radius)
            * velocity**3
            * (1 - gas.cp * wall / enthalpy)
        )

    def output(value: FloatArray) -> FloatResult:
        return return_float(value, scalar=scalar)

    return StagnationHeatingResult(
        model,
        output(flux),
        output(velocity),
        output(density),
        output(t0),
        output(ps),
        output(enthalpy),
        output(wall),
        gradient,
    )


def swept_cylinder_heat_flux(
    mach: ArrayLike,
    temperature: ArrayLike,
    pressure: ArrayLike,
    radius: ArrayLike,
    sweep_angle: ArrayLike,
    *,
    wall_temperature: ArrayLike,
    gas: PerfectGas = AIR,
) -> SweptCylinderHeatingResult:
    """Estimate laminar swept infinite-cylinder stagnation-line heat flux.

    Inputs and outputs use SI units, sweep is in radians, and heat flux is
    positive into the wall. Equation 41 is an air correlation; ``gas`` only
    converts the freestream state. The freestream must be supersonic.
    """
    (m, t, p, r, sweep, wall), scalar = _inputs(
        ("mach", mach),
        ("temperature", temperature),
        ("pressure", pressure),
        ("radius", radius),
        ("sweep_angle", sweep_angle),
        ("wall_temperature", wall_temperature),
    )
    if np.any((sweep < 0) | (sweep > np.pi / 2)):
        raise ValueError("sweep_angle must be between zero and pi/2")
    density, velocity = _freestream(m, t, p, gas)
    _tauber_warnings(velocity, wall, gas.cp * t)
    factor = 1 - TAUBER_SWEEP_FACTOR * np.sin(sweep) ** 2
    recovery = gas.cp * t + 0.5 * velocity**2 * factor
    flux = (
        TAUBER_CYLINDER_COEFFICIENT
        * np.sqrt(density / r)
        * factor
        * velocity**3
        * (1 - gas.cp * wall / recovery)
        * np.cos(sweep)
    )

    def output(value: FloatArray) -> FloatResult:
        return return_float(value, scalar=scalar)

    return SweptCylinderHeatingResult(
        output(flux),
        output(sweep),
        output(velocity),
        output(density),
        output(recovery),
        output(wall),
    )


def _attachment_line_edge(
    m: FloatArray,
    t: FloatArray,
    p: FloatArray,
    r: FloatArray,
    sweep: FloatArray,
    gas: PerfectGas,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return T_e, p_e, rho_e, W_e and circumferential gradient A."""
    _, velocity = _freestream(m, t, p, gas)
    edge_velocity = velocity * np.sin(sweep)
    normal_mach = m * np.cos(sweep)
    gamma = gas.heat_capacity_ratio
    factor = 1 + (gamma - 1) / 2 * normal_mach**2
    edge_temperature = t * factor
    edge_pressure = (
        p
        * factor ** (gamma / (gamma - 1))
        * np.asarray(normal_shock(normal_mach, gas).total_pressure_ratio)
    )
    edge_density = edge_pressure / (gas.specific_gas_constant * edge_temperature)
    gradient = np.sqrt(2 * (edge_pressure - p) / edge_density) / r
    return edge_temperature, edge_pressure, edge_density, edge_velocity, gradient


def _turbulent_attachment_line(
    edge_temperature: FloatArray,
    edge_pressure: FloatArray,
    edge_velocity: FloatArray,
    velocity_gradient: FloatArray,
    wall_temperature: FloatArray,
    *,
    prandtl_number: float,
    recovery_factor: float,
    gas: PerfectGas,
    viscosity_model: DynamicViscosityModel,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return q_w, T_aw, T*, Rbar*, St* and M_e from an edge state.

    Implements the Poll reference-temperature correlation as reproduced in
    NASA CR-181893 (5.9)--(5.10), using reference density at edge pressure.
    """
    recovery = edge_temperature + recovery_factor * edge_velocity**2 / (2 * gas.cp)
    reference = (
        edge_temperature
        + 0.1 * (wall_temperature - edge_temperature)
        + 0.6 * (recovery - edge_temperature)
    )
    density = edge_pressure / (gas.specific_gas_constant * reference)
    viscosity = viscosity_model.dynamic_viscosity(reference)
    kinematic_viscosity = viscosity / density
    reynolds = edge_velocity / np.sqrt(kinematic_viscosity * velocity_gradient)
    stanton = (
        POLL_TURBULENT_COEFFICIENT
        / (2 * prandtl_number ** (2 / 3))
        * reynolds**-POLL_TURBULENT_EXPONENT
    )
    flux = stanton * density * edge_velocity * gas.cp * (recovery - wall_temperature)
    edge_mach = edge_velocity / np.sqrt(
        gas.heat_capacity_ratio * gas.specific_gas_constant * edge_temperature
    )
    return flux, recovery, reference, reynolds, stanton, edge_mach


def swept_cylinder_turbulent_heat_flux(
    mach: ArrayLike,
    temperature: ArrayLike,
    pressure: ArrayLike,
    radius: ArrayLike,
    sweep_angle: ArrayLike,
    *,
    wall_temperature: ArrayLike,
    model: AttachmentLineHeatingModel = AttachmentLineHeatingModel.POLL,
    prandtl_number: float = 0.72,
    recovery_factor: float = DEFAULT_ATTACHMENT_LINE_RECOVERY_FACTOR,
    gas: PerfectGas = AIR,
    viscosity_model: DynamicViscosityModel = AIR_VISCOSITY,
) -> AttachmentLineHeatingResult:
    """Estimate fully turbulent heating at a long swept-cylinder attachment line.

    The user explicitly selects a smooth-wall, fully turbulent assumption;
    no transition or intermittency model is applied. Inputs use SI units,
    ``radius`` is the cylinder radius in m and ``sweep_angle`` is in radians
    (zero means the axis is perpendicular to the freestream). Require
    0 < sweep < pi/2 and mach*cos(sweep) > 1. Heat flux is positive into the
    wall and remains negative for a wall hotter than the recovery temperature.

    Modified-Newtonian edge properties stop only the normal velocity component.
    ``gas`` converts freestream and edge states, but the correlation is based
    on air data and has not been validated for arbitrary gases. This is the
    Poll reference-temperature correlation as reproduced in NASA CR-181893,
    not an independently checked transcription of Poll (1981).
    See :doc:`../models/aerodynamic_heating` for applicability limits.
    """
    if not isinstance(model, AttachmentLineHeatingModel):
        raise ValueError("model must be an AttachmentLineHeatingModel value")
    for name, value in (
        ("prandtl_number", prandtl_number),
        ("recovery_factor", recovery_factor),
    ):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and greater than zero")
    (m, t, p, r, sweep, wall), scalar = _inputs(
        ("mach", mach),
        ("temperature", temperature),
        ("pressure", pressure),
        ("radius", radius),
        ("sweep_angle", sweep_angle),
        ("wall_temperature", wall_temperature),
    )
    if np.any((sweep <= 0) | (sweep >= np.pi / 2)):
        raise ValueError(
            "sweep_angle must be greater than zero and less than pi/2 "
            "for the turbulent attachment line"
        )
    if np.any(m * np.cos(sweep) <= 1):
        raise ValueError(
            "normal Mach number mach*cos(sweep_angle) must be greater than one"
        )
    density, velocity = _freestream(m, t, p, gas)
    te, pe, rho_e, we, gradient = _attachment_line_edge(m, t, p, r, sweep, gas)
    flux, recovery, reference, reynolds, stanton, edge_mach = (
        _turbulent_attachment_line(
            te,
            pe,
            we,
            gradient,
            wall,
            prandtl_number=prandtl_number,
            recovery_factor=recovery_factor,
            gas=gas,
            viscosity_model=viscosity_model,
        )
    )
    total_temperature = t * (1 + (gas.heat_capacity_ratio - 1) / 2 * m**2)
    if np.any(total_temperature > PERFECT_GAS_STAGNATION_TEMPERATURE_LIMIT):
        warnings.warn(
            "Turbulent attachment-line edge properties use a perfect gas above "
            "2000 K total temperature; dissociation is not modelled",
            ApplicabilityWarning,
            stacklevel=2,
        )
    if np.any(
        (edge_mach < BENARD_MINIMUM_EDGE_MACH) | (edge_mach > BENARD_MAXIMUM_EDGE_MACH)
    ):
        warnings.warn(
            "Turbulent attachment-line edge Mach number is outside "
            "0.35 <= M_e <= 6.89, the turbulent heat-transfer database range "
            "of Benard et al. (2006)",
            ApplicabilityWarning,
            stacklevel=2,
        )
    wall_ratio = wall / total_temperature
    if np.any(
        (wall_ratio < BENARD_MINIMUM_WALL_TOTAL_RATIO)
        | (wall_ratio > BENARD_MAXIMUM_WALL_TOTAL_RATIO)
    ):
        warnings.warn(
            "Turbulent attachment-line wall-to-total temperature ratio is outside "
            "0.14 <= T_w/T_0 <= 0.96, the turbulent heat-transfer database range "
            "of Benard et al. (2006)",
            ApplicabilityWarning,
            stacklevel=2,
        )

    def output(value: FloatArray) -> FloatResult:
        return return_float(value, scalar=scalar)

    return AttachmentLineHeatingResult(
        model=model,
        heat_flux=output(flux),
        sweep_angle=output(sweep),
        freestream_velocity=output(velocity),
        freestream_density=output(density),
        edge_temperature=output(te),
        edge_pressure=output(pe),
        edge_density=output(rho_e),
        edge_velocity=output(we),
        edge_mach=output(edge_mach),
        velocity_gradient=output(gradient),
        recovery_temperature=output(recovery),
        reference_temperature=output(reference),
        reynolds_number=output(reynolds),
        stanton_number=output(stanton),
        wall_temperature=output(wall),
    )


__all__ = [
    "DEFAULT_ATTACHMENT_LINE_RECOVERY_FACTOR",
    "POLL_TURBULENT_COEFFICIENT",
    "POLL_TURBULENT_EXPONENT",
    "TAUBER_CYLINDER_COEFFICIENT",
    "TAUBER_SPHERE_COEFFICIENT",
    "TAUBER_SWEEP_FACTOR",
    "AttachmentLineHeatingModel",
    "AttachmentLineHeatingResult",
    "StagnationHeatingModel",
    "StagnationHeatingResult",
    "SweptCylinderHeatingResult",
    "sphere_stagnation_heat_flux",
    "swept_cylinder_heat_flux",
    "swept_cylinder_turbulent_heat_flux",
]
