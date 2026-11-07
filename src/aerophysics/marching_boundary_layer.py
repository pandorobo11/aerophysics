"""Boundary layers with streamwise-varying edge conditions.

The solver marches the thin-layer equations from a sharp leading edge along the
wall, so each station inherits the upstream pressure and wall-temperature
history. The edge state is an input: this module does not compute the outer
flow or represent a body shape.

All inputs and outputs use SI units. ``distance`` is measured along the wall
from the leading edge, the skin-friction coefficient uses the local edge dynamic
pressure, and heat flux into the wall is positive.

The current scope is planar laminar or SA-noft2 flow of a perfect gas with an
isentropic edge. Separation stops the march because the downstream march is not
valid past it.
"""

import warnings
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import PchipInterpolator

from aerophysics._array import FloatArray, as_float_array
from aerophysics._boundary_layer_march import (
    Discretization,
    MarchSolution,
    SeparationError,
    StationEdge,
    eta_grid,
    march,
)
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics._native_march import native_march
from aerophysics.boundary_layer import BoundaryLayerRegime
from aerophysics.exceptions import (
    ApplicabilityWarning,
    BoundaryLayerSeparationError,
    ModelRangeError,
)
from aerophysics.gas import AIR, PerfectGas
from aerophysics.transport import AIR_VISCOSITY, DynamicViscosityModel

_LAMINAR_ETA_FIRST = 1e-3
_LAMINAR_ETA_RATIO = 1.015
_LAMINAR_ETA_SPACING = 0.05
_LAMINAR_ETA_MAXIMUM = 16.0

_SA_ETA_FIRST = 7e-4
_SA_ETA_RATIO = 1.04
_SA_ETA_SPACING = 1.0
_SA_ETA_MINIMUM_EXTENT = 150.0
_SA_ETA_REFERENCE_REYNOLDS = 1e8
_SA_ETA_REYNOLDS_EXPONENT = 0.42


def _turbulent_eta_grid(reynolds_equivalent: float) -> FloatArray:
    """Size the outer domain before marching, preserving the low-Re grid."""
    scale = max(1.0, reynolds_equivalent / _SA_ETA_REFERENCE_REYNOLDS) ** (
        _SA_ETA_REYNOLDS_EXPONENT
    )
    # Do not change the reference grid for roundoff at Re_equivalent = 1e8.
    if scale <= 1.0 + 32.0 * np.finfo(np.float64).eps:
        scale = 1.0
    return eta_grid(
        _SA_ETA_FIRST,
        _SA_ETA_RATIO,
        _SA_ETA_SPACING * scale,
        _SA_ETA_MINIMUM_EXTENT * scale,
    )


_PANELS_PER_INTERVAL = 16
_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(5)

_EDGE_CONSISTENCY_TOLERANCE = 1e-6
_VISCOUS_INTERACTION_THRESHOLD = 1.0
_PRESSURE_LENGTH_THRESHOLD = 1.0
_VELOCITY_THRESHOLD = 0.99


def _read_only(values: FloatArray) -> FloatArray:
    result = np.array(values, dtype=np.float64, copy=True)
    result.flags.writeable = False
    return result


def _positive_scalar(value: float, *, name: str) -> float:
    array, scalar = as_float_array(value, name=name)
    if not scalar or float(array) <= 0.0:
        raise ValueError(f"{name} must be a scalar greater than zero")
    return float(array)


def _distance_array(distance: ArrayLike) -> FloatArray:
    array, scalar = as_float_array(distance, name="distance")
    if scalar or array.ndim != 1 or array.size < 2:
        raise ValueError("distance must be a one-dimensional array of two or more")
    if array[0] != 0.0:
        raise ValueError("distance must start at the leading edge, 0")
    if np.any(np.diff(array) <= 0.0):
        raise ValueError("distance must increase strictly")
    return np.array(array, dtype=np.float64, copy=True)


def _station_array(value: ArrayLike, count: int, *, name: str) -> FloatArray:
    array, scalar = as_float_array(value, name=name)
    if scalar:
        array = np.full(count, float(array), dtype=np.float64)
    if array.shape != (count,):
        raise ValueError(f"{name} must match the shape of distance")
    if np.any(array <= 0.0):
        raise ValueError(f"{name} must be greater than zero")
    return np.array(array, dtype=np.float64, copy=True)


def _isentropic_state(
    pressure: FloatArray,
    total_pressure: float,
    total_temperature: float,
    gas: PerfectGas,
) -> tuple[FloatArray, FloatArray]:
    gamma = gas.heat_capacity_ratio
    temperature = total_temperature * (pressure / total_pressure) ** (
        (gamma - 1.0) / gamma
    )
    velocity = np.sqrt(2.0 * gas.cp * (total_temperature - temperature))
    return temperature, velocity


@dataclass(frozen=True, slots=True)
class BoundaryLayerEdge:
    """Isentropic edge state along the wall, starting at the leading edge.

    The governing equations require ``dp = -rho U dU`` and a constant total
    enthalpy along the edge, so the three distributions are not independent.
    Prefer :meth:`from_pressure`. Direct construction checks these relations.

    Parameters
    ----------
    distance:
        Distance along the wall from the leading edge in m. It starts at zero
        and increases strictly.
    pressure:
        Edge static pressure in Pa.
    temperature:
        Edge static temperature in K.
    velocity:
        Edge velocity in m/s, greater than zero at every station.
    gas:
        Calorically perfect gas.
    viscosity_model:
        Dynamic-viscosity model.

    Raises
    ------
    ValueError
        If an array is invalid or the edge state is not isentropic.
    """

    distance: FloatArray
    pressure: FloatArray
    temperature: FloatArray
    velocity: FloatArray
    gas: PerfectGas = AIR
    viscosity_model: DynamicViscosityModel = AIR_VISCOSITY

    def __post_init__(self) -> None:
        distance = _distance_array(self.distance)
        count = distance.size
        pressure = _station_array(self.pressure, count, name="pressure")
        temperature = _station_array(self.temperature, count, name="temperature")
        velocity = _station_array(self.velocity, count, name="velocity")
        total_temperature = temperature[0] + velocity[0] ** 2 / (2.0 * self.gas.cp)
        gamma = self.gas.heat_capacity_ratio
        total_pressure = pressure[0] * (total_temperature / temperature[0]) ** (
            gamma / (gamma - 1.0)
        )
        if np.any(pressure >= total_pressure):
            raise ValueError(
                "pressure must stay below the total pressure of the first station; "
                "the edge state must be isentropic with constant total enthalpy"
            )
        expected_temperature, expected_velocity = _isentropic_state(
            pressure, total_pressure, total_temperature, self.gas
        )
        if np.any(
            np.abs(temperature / expected_temperature - 1.0)
            > _EDGE_CONSISTENCY_TOLERANCE
        ) or np.any(
            np.abs(velocity / expected_velocity - 1.0) > _EDGE_CONSISTENCY_TOLERANCE
        ):
            raise ValueError(
                "the edge state must be isentropic with constant total enthalpy; "
                "use BoundaryLayerEdge.from_pressure"
            )
        for name, values in (
            ("distance", distance),
            ("pressure", pressure),
            ("temperature", temperature),
            ("velocity", velocity),
        ):
            object.__setattr__(self, name, _read_only(values))

    @classmethod
    def from_pressure(
        cls,
        distance: ArrayLike,
        pressure: ArrayLike,
        *,
        total_pressure: float,
        total_temperature: float,
        gas: PerfectGas = AIR,
        viscosity_model: DynamicViscosityModel = AIR_VISCOSITY,
    ) -> "BoundaryLayerEdge":
        """Build the edge state from a pressure distribution.

        Parameters
        ----------
        distance:
            Distance along the wall from the leading edge in m.
        pressure:
            Edge static pressure in Pa, below ``total_pressure``.
        total_pressure:
            Edge total pressure in Pa. Downstream of a shock this is the
            post-shock value.
        total_temperature:
            Edge total temperature in K.
        gas:
            Calorically perfect gas.
        viscosity_model:
            Dynamic-viscosity model.

        Returns
        -------
        BoundaryLayerEdge
            Edge state whose temperature and velocity follow from isentropic
            relations.

        Raises
        ------
        ValueError
            If an input is invalid or a pressure reaches the total pressure.
        """
        stations = _distance_array(distance)
        pressures = _station_array(pressure, stations.size, name="pressure")
        total_p = _positive_scalar(total_pressure, name="total_pressure")
        total_t = _positive_scalar(total_temperature, name="total_temperature")
        if np.any(pressures >= total_p):
            raise ValueError("pressure must stay below total_pressure")
        temperature, velocity = _isentropic_state(pressures, total_p, total_t, gas)
        return cls(stations, pressures, temperature, velocity, gas, viscosity_model)


@dataclass(frozen=True, slots=True)
class MarchingBoundaryLayerResult:
    """Boundary layer at each edge station downstream of the leading edge.

    The leading edge itself is excluded because the wall shear and heat flux are
    singular there. One-dimensional arrays have the shape ``(N_s,)`` and
    two-dimensional arrays the shape ``(N_s, N_eta)``.

    ``displacement_thickness`` and ``momentum_thickness`` integrate to the outer
    numerical boundary. ``wall_heat_flux`` follows from the computed temperature
    gradient at the wall, not from a Reynolds analogy; it is zero for an
    adiabatic wall. ``pressure_gradient_parameter`` is
    ``beta = 2 xi d(ln U_e)/d(xi)``, positive for accelerating flow.
    ``nu_tilde_ratio`` is the read-only SA working variable divided by local
    edge kinematic viscosity, or None for laminar flow.
    Its outer profile can oscillate and is not grid-converged in strong
    expansions or with a small ``edge_sa_ratio``. Wall quantities and
    ``delta99`` were insensitive to outer-grid refinement in the measured cases.
    """

    distance: FloatArray
    reynolds_number: FloatArray
    boundary_layer_thickness: FloatArray
    displacement_thickness: FloatArray
    momentum_thickness: FloatArray
    local_skin_friction_coefficient: FloatArray
    wall_shear_stress: FloatArray
    wall_temperature: FloatArray
    wall_heat_flux: FloatArray
    pressure_gradient_parameter: FloatArray
    eta: FloatArray
    wall_distance: FloatArray
    velocity_ratio: FloatArray
    temperature_ratio: FloatArray
    nu_tilde_ratio: FloatArray | None = None


class _EdgeInterpolant:
    """Monotone interpolation of the edge state and of the ``xi`` integral."""

    def __init__(
        self,
        edge: BoundaryLayerEdge,
        wall_temperature: FloatArray | None,
        edge_sa_ratio: float | None = None,
    ) -> None:
        self.edge_sa_ratio = edge_sa_ratio
        self.initial_nu = float(
            edge.viscosity_model.dynamic_viscosity(edge.temperature[0])
        ) / (edge.pressure[0] / (edge.gas.specific_gas_constant * edge.temperature[0]))
        self.gas = edge.gas
        self.viscosity_model = edge.viscosity_model
        self.total_temperature = edge.temperature[0] + edge.velocity[0] ** 2 / (
            2.0 * edge.gas.cp
        )
        gamma = edge.gas.heat_capacity_ratio
        self.total_pressure = edge.pressure[0] * (
            self.total_temperature / edge.temperature[0]
        ) ** (gamma / (gamma - 1.0))
        self.log_pressure = PchipInterpolator(edge.distance, np.log(edge.pressure))
        self.log_pressure_slope = self.log_pressure.derivative()
        self.wall_temperature = (
            None
            if wall_temperature is None
            else PchipInterpolator(edge.distance, wall_temperature)
        )
        panels = np.concatenate(
            [
                np.linspace(left, right, _PANELS_PER_INTERVAL, endpoint=False)
                for left, right in zip(
                    edge.distance[:-1], edge.distance[1:], strict=True
                )
            ]
            + [edge.distance[-1:]]
        )
        half = (panels[1:] - panels[:-1]) / 2.0
        nodes = panels[:-1, None] + half[:, None] * (_GAUSS_NODES + 1.0)
        flux = self._flux(nodes)
        # Copy each row to match the scalar integral's aligned flux allocation.
        # Unaligned five-element views can change np.dot's summation order.
        increments = np.array(
            [
                float(h * np.dot(_GAUSS_WEIGHTS, row.copy()))
                for h, row in zip(half, flux, strict=True)
            ]
        )
        self.panels = panels
        self.xi_at_panels = np.concatenate(([0.0], np.cumsum(increments)))

    def _flux(self, distance: FloatArray) -> FloatArray:
        """Return ``rho_e mu_e U_e``."""
        pressure = np.exp(self.log_pressure(distance))
        temperature, velocity = _isentropic_state(
            pressure, self.total_pressure, self.total_temperature, self.gas
        )
        density = pressure / (self.gas.specific_gas_constant * temperature)
        viscosity = np.asarray(self.viscosity_model.dynamic_viscosity(temperature))
        return np.asarray(density * viscosity * velocity, dtype=np.float64)

    def _integral(self, left: float, right: float) -> float:
        half = (right - left) / 2.0
        nodes = left + half * (_GAUSS_NODES + 1.0)
        return float(half * np.dot(_GAUSS_WEIGHTS, self._flux(nodes)))

    def xi(self, distance: float) -> float:
        """Return ``int_0^s rho_e mu_e U_e ds``."""
        index = int(np.searchsorted(self.panels, distance, side="right")) - 1
        index = min(max(index, 0), len(self.panels) - 1)
        left = float(self.panels[index])
        increment = self._integral(left, distance) if distance > left else 0.0
        return float(self.xi_at_panels[index]) + increment

    def state(self, distance: float) -> StationEdge:
        """Return the transformed-equation parameters at one position."""
        gas = self.gas
        pressure = float(np.exp(self.log_pressure(distance)))
        temperature = self.total_temperature * (pressure / self.total_pressure) ** (
            (gas.heat_capacity_ratio - 1.0) / gas.heat_capacity_ratio
        )
        velocity = float(np.sqrt(2.0 * gas.cp * (self.total_temperature - temperature)))
        density = pressure / (gas.specific_gas_constant * temperature)
        viscosity = float(self.viscosity_model.dynamic_viscosity(temperature))
        xi = self.xi(distance)
        # dp = -rho U dU on the edge streamline.
        log_velocity_slope = (
            -pressure
            * float(self.log_pressure_slope(distance))
            / (density * velocity**2)
        )
        beta = 2.0 * xi * log_velocity_slope / (density * viscosity * velocity)
        heating = velocity**2 / (gas.cp * temperature)
        wall_ratio = (
            None
            if self.wall_temperature is None
            else float(self.wall_temperature(distance)) / temperature
        )
        if self.edge_sa_ratio is None:
            growth = 0.0
            sa_ratio = 0.0
        else:
            # Only SA needs this derivative, including at zero edge SA ratio.
            # The viscosity model must accept T_e * exp(+/-1e-5).
            delta = 1e-5
            mu_plus = float(
                self.viscosity_model.dynamic_viscosity(temperature * np.exp(delta))
            )
            mu_minus = float(
                self.viscosity_model.dynamic_viscosity(temperature * np.exp(-delta))
            )
            exponent = np.log(mu_plus / mu_minus) / (2.0 * delta)
            gamma = gas.heat_capacity_ratio
            log_nu_slope = (
                ((gamma - 1.0) * exponent - 1.0)
                / gamma
                * float(self.log_pressure_slope(distance))
            )
            growth = 2.0 * xi * log_nu_slope / (density * viscosity * velocity)
            sa_ratio = self.edge_sa_ratio * self.initial_nu * density / viscosity
        return StationEdge(
            xi,
            velocity,
            beta,
            heating,
            temperature,
            wall_ratio,
            float(np.sqrt(2.0 * xi) / viscosity),
            float(growth),
            sa_ratio,
        )


def _first_crossing(
    velocity: FloatArray, coordinate: FloatArray, level: float
) -> float:
    """Interpolate the coordinate where ``velocity`` first reaches ``level``."""
    above = int(np.argmax(velocity >= level))
    fraction = (level - velocity[above - 1]) / (velocity[above] - velocity[above - 1])
    return float(
        coordinate[above - 1] + fraction * (coordinate[above] - coordinate[above - 1])
    )


def _result(
    edge: BoundaryLayerEdge,
    interpolant: _EdgeInterpolant,
    solution: MarchSolution,
    prandtl_number: float,
    edge_values: FloatArray | None = None,
) -> MarchingBoundaryLayerResult:
    gas = edge.gas
    count = solution.distance.size
    eta = solution.eta
    velocity_ratio = solution.velocity_ratio
    temperature_ratio = solution.temperature_ratio
    spacing = np.diff(eta)
    stretched = np.zeros_like(temperature_ratio)
    stretched[:, 1:] = np.cumsum(
        spacing * (temperature_ratio[:, 1:] + temperature_ratio[:, :-1]) / 2.0, axis=1
    )

    def field() -> FloatArray:
        return np.empty(count, dtype=np.float64)

    reynolds, thickness, displacement, momentum = field(), field(), field(), field()
    friction, shear, wall_temperature = field(), field(), field()
    heat_flux, beta = field(), field()
    wall_distance = np.empty_like(temperature_ratio)
    for index, distance in enumerate(float(value) for value in solution.distance):
        if edge_values is None:
            state = interpolant.state(distance)
            pressure = float(np.exp(interpolant.log_pressure(distance)))
        else:
            xi, velocity, beta_value, temperature, pressure = edge_values[index, :5]
            state = StationEdge(xi, velocity, beta_value, 0.0, temperature, None)
        temperature, xi = state.temperature, state.xi
        density = pressure / (gas.specific_gas_constant * temperature)
        viscosity = float(edge.viscosity_model.dynamic_viscosity(temperature))
        root = np.sqrt(2.0 * xi)
        scale = root / (state.velocity * density)
        u, g = velocity_ratio[index], temperature_ratio[index]
        wall = g[0] * temperature
        wall_viscosity = float(edge.viscosity_model.dynamic_viscosity(wall))
        wall_c = wall_viscosity / viscosity / g[0]
        reynolds[index] = density * state.velocity * distance / viscosity
        wall_distance[index] = scale * stretched[index]
        thickness[index] = _first_crossing(u, wall_distance[index], _VELOCITY_THRESHOLD)
        displacement[index] = scale * float(np.trapezoid(g - u, eta))
        momentum[index] = scale * float(np.trapezoid(u * (1.0 - u), eta))
        friction[index] = (
            2.0 * wall_c * viscosity * solution.wall_velocity_slope[index] / root
        )
        shear[index] = friction[index] * 0.5 * density * state.velocity**2
        wall_temperature[index] = wall
        heat_flux[index] = (
            0.0
            if interpolant.wall_temperature is None
            else wall_c
            * density
            * viscosity
            * state.velocity
            * gas.cp
            * temperature
            * solution.wall_temperature_slope[index]
            / (prandtl_number * root)
        )
        beta[index] = state.beta
    return MarchingBoundaryLayerResult(
        distance=_read_only(solution.distance),
        reynolds_number=_read_only(reynolds),
        boundary_layer_thickness=_read_only(thickness),
        displacement_thickness=_read_only(displacement),
        momentum_thickness=_read_only(momentum),
        local_skin_friction_coefficient=_read_only(friction),
        wall_shear_stress=_read_only(shear),
        wall_temperature=_read_only(wall_temperature),
        wall_heat_flux=_read_only(heat_flux),
        pressure_gradient_parameter=_read_only(beta),
        eta=_read_only(eta),
        wall_distance=_read_only(wall_distance),
        velocity_ratio=_read_only(velocity_ratio),
        temperature_ratio=_read_only(temperature_ratio),
        nu_tilde_ratio=None
        if solution.nu_tilde_ratio is None
        else _read_only(solution.nu_tilde_ratio),
    )


def _warn_applicability(
    edge: BoundaryLayerEdge,
    interpolant: _EdgeInterpolant,
    result: MarchingBoundaryLayerResult,
) -> None:
    if result.distance.size == 0:
        return
    pressure_length = np.abs(interpolant.log_pressure_slope(result.distance)) * (
        result.boundary_layer_thickness
    )
    if np.max(pressure_length) > _PRESSURE_LENGTH_THRESHOLD:
        warnings.warn(
            "the edge pressure changes over a distance comparable to the boundary "
            "layer thickness, which violates the thin-layer approximation (maximum "
            f"delta99 |d ln p/ds| = {np.max(pressure_length):.6g})",
            ApplicabilityWarning,
            stacklevel=3,
        )
    pressure = np.exp(interpolant.log_pressure(result.distance))
    gas = edge.gas
    temperature, velocity = _isentropic_state(
        pressure, interpolant.total_pressure, interpolant.total_temperature, gas
    )
    mach = velocity / np.sqrt(
        gas.heat_capacity_ratio * gas.specific_gas_constant * temperature
    )
    viscosity = np.asarray(edge.viscosity_model.dynamic_viscosity(temperature))
    wall_viscosity = np.asarray(
        edge.viscosity_model.dynamic_viscosity(result.wall_temperature)
    )
    wall_c = temperature / result.wall_temperature * wall_viscosity / viscosity
    if result.nu_tilde_ratio is not None:
        wall_density = pressure / (gas.specific_gas_constant * result.wall_temperature)
        y_plus = (
            result.wall_distance[:, 1]
            * np.sqrt(wall_density * np.maximum(result.wall_shear_stress, 0.0))
            / wall_viscosity
        )
        if np.max(y_plus) > 1.0:
            warnings.warn(
                "the first wall-normal node has y+ > 1 "
                f"(maximum y+ = {np.max(y_plus):.6g})",
                ApplicabilityWarning,
                stacklevel=3,
            )
    interaction = mach**3 * np.sqrt(wall_c / result.reynolds_number)
    if np.max(interaction) > _VISCOUS_INTERACTION_THRESHOLD:
        warnings.warn(
            "laminar viscous interaction: leading-edge induced pressure is not "
            "represented by the supplied edge state "
            f"(maximum chi_bar = {np.max(interaction):.6g})",
            ApplicabilityWarning,
            stacklevel=3,
        )


def marching_boundary_layer(
    edge: BoundaryLayerEdge,
    *,
    regime: BoundaryLayerRegime,
    wall_temperature: ArrayLike | None = None,
    prandtl_number: float = 0.72,
    turbulent_prandtl_number: float = 0.9,
    edge_sa_ratio: float = 3.0,
) -> MarchingBoundaryLayerResult:
    """March a boundary layer along a wall from a sharp leading edge.

    Parameters
    ----------
    edge:
        Edge state along the wall. Resolve pressure changes smoothly, because
        the march follows a monotone interpolation of the supplied pressure.
    regime:
        ``LAMINAR`` or ``TURBULENT`` (SA-noft2). ``TRANSITIONAL`` is not
        implemented. SA-noft2 develops its own laminar start and transition
        near the leading edge; this is not a transition prediction.
    wall_temperature:
        Wall temperature in K, as a scalar or an array with the shape of
        ``edge.distance``. None selects an adiabatic wall.
    prandtl_number:
        Constant Prandtl number.
    turbulent_prandtl_number:
        Constant turbulent Prandtl number, greater than zero. Validated in
        every regime and used only for ``TURBULENT``.
    edge_sa_ratio:
        Initial ``nu_tilde/nu_e`` at the edge, at least zero. Dimensional
        ``nu_tilde`` is conserved along the edge streamline. Zero gives the
        laminar limit. Validated in every regime and used only for ``TURBULENT``.

    Returns
    -------
    MarchingBoundaryLayerResult
        Boundary layer at the edge stations downstream of the leading edge.

    Raises
    ------
    ValueError
        If an input is invalid.
    ModelRangeError
        If the requested regime is not implemented.
    BoundaryLayerSeparationError
        If the wall shear collapses in an adverse pressure gradient. The
        exception holds the separation distance and the upstream result.
    BoundaryLayerConvergenceError
        If the march cannot resolve the solution.

    Warns
    -----
    ApplicabilityWarning
        If the pressure changes over a distance comparable to the boundary-layer
        thickness, laminar viscous interaction is strong, or the turbulent
        first wall-normal node has ``y+ > 1`` using local wall properties.

    Notes
    -----
    Turbulent flow requires the viscosity model to accept ``T_e * exp(+/-1e-5)``
    for an edge viscosity derivative, even when ``edge_sa_ratio=0``. Near a
    model's temperature-range limit, these evaluations can raise its range
    exception or warning. Laminar flow does not evaluate this SA derivative.

    Examples
    --------
    >>> import numpy as np
    >>> from aerophysics import (
    ...     BoundaryLayerEdge,
    ...     BoundaryLayerRegime,
    ...     marching_boundary_layer,
    ... )
    >>> edge = BoundaryLayerEdge.from_pressure(
    ...     [0.0, 0.5, 1.0],
    ...     np.full(3, 20000.0),
    ...     total_pressure=156000.0,
    ...     total_temperature=500.0,
    ... )
    >>> layer = marching_boundary_layer(edge, regime=BoundaryLayerRegime.LAMINAR)
    >>> layer.distance.tolist()
    [0.5, 1.0]
    """
    turbulent_prandtl = _positive_scalar(
        turbulent_prandtl_number, name="turbulent_prandtl_number"
    )
    sa_array, sa_scalar = as_float_array(edge_sa_ratio, name="edge_sa_ratio")
    if not sa_scalar or float(sa_array) < 0.0:
        raise ValueError("edge_sa_ratio must be a scalar at least zero")
    sa_ratio = float(sa_array)
    if (
        regime is not BoundaryLayerRegime.LAMINAR
        and regime is not BoundaryLayerRegime.TURBULENT
    ):
        raise ModelRangeError(
            "marching_boundary_layer implements BoundaryLayerRegime.LAMINAR "
            "and TURBULENT"
        )
    prandtl = _positive_scalar(prandtl_number, name="prandtl_number")
    wall = (
        None
        if wall_temperature is None
        else _station_array(
            wall_temperature, edge.distance.size, name="wall_temperature"
        )
    )
    interpolant = _EdgeInterpolant(
        edge, wall, sa_ratio if regime is BoundaryLayerRegime.TURBULENT else None
    )
    eta = eta_grid(
        _LAMINAR_ETA_FIRST,
        _LAMINAR_ETA_RATIO,
        _LAMINAR_ETA_SPACING,
        _LAMINAR_ETA_MAXIMUM,
    )
    discretization = (
        TurbulentDiscretization(
            _turbulent_eta_grid(
                max(
                    interpolant.state(float(s)).reynolds_root for s in edge.distance[1:]
                )
                ** 2
                / 2.0
            ),
            prandtl,
            turbulent_prandtl,
            edge.viscosity_model,
        )
        if regime is BoundaryLayerRegime.TURBULENT
        else Discretization(eta, prandtl, edge.viscosity_model)
    )
    edge_values = None
    try:
        native = (
            native_march(
                interpolant, edge.distance[1:], discretization, edge.viscosity_model
            )
            if regime is BoundaryLayerRegime.LAMINAR
            else None
        )
        if native is None:
            solution = march(
                interpolant.state,
                edge.distance[1:],
                discretization,
            )
        else:
            solution, edge_values = native
    except SeparationError as separation:
        upstream = _result(edge, interpolant, separation.solution, prandtl)
        raise BoundaryLayerSeparationError(
            "the wall shear collapses at distance "
            f"{separation.distance:.6g} m in an adverse pressure gradient",
            separation_distance=separation.distance,
            result=upstream,
        ) from separation
    result = _result(edge, interpolant, solution, prandtl, edge_values)
    _warn_applicability(edge, interpolant, result)
    return result


__all__ = [
    "BoundaryLayerEdge",
    "MarchingBoundaryLayerResult",
    "marching_boundary_layer",
]
