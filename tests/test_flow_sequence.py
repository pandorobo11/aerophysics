"""Conservation, independent closed forms, handoffs and failure boundaries."""

import math
from typing import Any, cast

import pytest

from aerophysics.exceptions import ExpansionConvergenceError, ModelRangeError
from aerophysics.flow_sequence import (
    ConicalShockStep,
    ExpansionStep,
    FlowProperties,
    FlowSequenceResult,
    FlowState,
    FlowStep,
    IsentropicStep,
    NormalShockStep,
    ObliqueShockStep,
    flow_properties,
    solve_flow_sequence,
)
from aerophysics.gas import AIR
from aerophysics.isentropic import MachBranch
from aerophysics.real_gas import AIR_HARMONIC_OSCILLATOR
from aerophysics.shocks import ShockBranch, ShockGasModel, conical_shock, oblique_shock
from aerophysics.thermochemistry import AIR_NASA7, AIR_NASA9


def exit_properties(result: FlowSequenceResult, index: int) -> FlowProperties:
    properties = result.steps[index].properties
    assert properties is not None
    return properties


@pytest.mark.parametrize("gas", [AIR, AIR_NASA7, AIR_NASA9, AIR_HARMONIC_OSCILLATOR])
def test_sequence_conservation_and_isentropic_equivalence(gas: ShockGasModel) -> None:
    initial = FlowState(23000, 500, 2)
    result = solve_flow_sequence(
        initial,
        [ObliqueShockStep(math.radians(5)), ExpansionStep(math.radians(3))],
        gas,
    )
    assert result.complete
    previous = result.initial
    for stage in result.steps:
        current = stage.properties
        assert current is not None
        assert current.total_pressure is not None
        assert previous.total_pressure is not None
        assert result.initial.total_pressure is not None
        assert stage.entropy_change is not None
        assert current.total_enthalpy == pytest.approx(
            previous.total_enthalpy, rel=2e-8
        )
        assert current.total_pressure <= previous.total_pressure * (1 + 1e-9)
        assert stage.total_pressure_recovery == pytest.approx(
            current.total_pressure / result.initial.total_pressure, rel=2e-8
        )
        assert stage.entropy_change >= 0
        previous = current
    assert exit_properties(result, 1).total_pressure == pytest.approx(
        exit_properties(result, 0).total_pressure, rel=2e-8
    )
    target = exit_properties(
        solve_flow_sequence(initial, [IsentropicStep(2.2)], gas), 0
    )
    for step in [
        IsentropicStep(target.state.pressure / initial.pressure, "pressure_ratio"),
        IsentropicStep(
            result.initial.mass_flux / target.mass_flux,
            "area_ratio",
            MachBranch.SUPERSONIC,
        ),
    ]:
        replay = solve_flow_sequence(initial, [step], gas)
        assert replay.complete
        assert exit_properties(replay, 0).state.mach == pytest.approx(2.2, rel=2e-9)
        assert exit_properties(replay, 0).state.temperature == pytest.approx(
            target.state.temperature, rel=2e-9
        )


def test_normal_shock_closed_form_and_mass_momentum() -> None:
    result = solve_flow_sequence(FlowState(10000, 300, 3), [NormalShockStep()])
    before, after = result.initial, exit_properties(result, 0)
    assert after.state.pressure / before.state.pressure == pytest.approx(31 / 3)
    assert after.density / before.density == pytest.approx(27 / 7)
    assert after.state.mach == pytest.approx(math.sqrt(7 / 31))
    assert after.mass_flux == pytest.approx(before.mass_flux)
    assert after.state.pressure + after.density * after.velocity**2 == pytest.approx(
        before.state.pressure + before.density * before.velocity**2
    )
    assert after.total_enthalpy == pytest.approx(before.total_enthalpy)


@pytest.mark.parametrize("gas", [AIR, AIR_NASA7])
def test_cone_surface_handoff_and_strong_shock(gas: ShockGasModel) -> None:
    initial = FlowState(10000, 300, 3)
    cone = conical_shock(3, math.radians(10), gas, upstream_temperature=300)
    result = solve_flow_sequence(
        initial, [ConicalShockStep(math.radians(10)), NormalShockStep()], gas
    )
    assert result.complete
    assert exit_properties(result, 0).state.mach == pytest.approx(cone.surface_mach)
    assert result.steps[0].post_shock_mach == pytest.approx(cone.post_shock_mach)
    assert exit_properties(result, 1).total_enthalpy == pytest.approx(
        result.initial.total_enthalpy, rel=2e-8
    )
    strong = solve_flow_sequence(
        initial, [ObliqueShockStep(0.1, ShockBranch.STRONG)], gas
    )
    direct = oblique_shock(3, 0.1, ShockBranch.STRONG, gas, upstream_temperature=300)
    assert exit_properties(strong, 0).state.mach == pytest.approx(
        direct.downstream_mach
    )


def test_initialization_and_empty() -> None:
    initial = FlowState.from_atmosphere(0, 2)
    assert initial.pressure == pytest.approx(101325)
    assert initial.temperature == pytest.approx(288.15)
    props = flow_properties(initial)
    assert (
        FlowState.from_velocity(initial.pressure, initial.temperature, props.velocity)
        == initial
    )
    assert props.total_pressure is not None
    assert props.total_temperature is not None
    restored = FlowState.from_total(props.total_pressure, props.total_temperature, 2)
    assert restored.pressure == pytest.approx(initial.pressure)
    assert restored.temperature == pytest.approx(initial.temperature)
    assert solve_flow_sequence(initial, []).complete
    assert flow_properties(
        FlowState(10000, 500, 0), AIR_NASA7
    ).total_temperature == pytest.approx(500)


@pytest.mark.parametrize(
    "state", [(0, 300, 2), (1, -1, 2), (1, 300, -1), (1, 300, math.nan)]
)
def test_invalid_initial(state: tuple[float, float, float]) -> None:
    with pytest.raises(ValueError):
        FlowState(*state)


@pytest.mark.parametrize("speed", [-1, math.nan])
def test_invalid_velocity(speed: float) -> None:
    with pytest.raises(ValueError):
        FlowState.from_velocity(1, 300, speed)


def test_invalid_temperature_and_range() -> None:
    with pytest.raises(ValueError):
        FlowState.from_velocity(1, -1, 3)
    with pytest.raises(ModelRangeError):
        solve_flow_sequence(FlowState(10000, 100, 2), [], AIR_NASA7)
    with pytest.raises(ModelRangeError):
        FlowState.from_atmosphere(90000, 2)


@pytest.mark.parametrize(
    "step,status",
    [
        (ObliqueShockStep(1), "detached"),
        (ExpansionStep(-1), "invalid_input"),
        (ObliqueShockStep(math.nan), "invalid_input"),
        (IsentropicStep(-1), "invalid_input"),
        (IsentropicStep(0, "pressure_ratio"), "invalid_input"),
        (IsentropicStep(0.01, "area_ratio", MachBranch.SUBSONIC), "invalid_input"),
        (IsentropicStep(1, "area_ratio"), "invalid_input"),
        (IsentropicStep(1, cast(Any, "unknown")), "invalid_input"),
    ],
)
def test_failures_preserve_prefix(step: FlowStep, status: str) -> None:
    result = solve_flow_sequence(
        FlowState(10000, 300, 3), [ExpansionStep(0), step, NormalShockStep()]
    )
    assert not result.complete
    assert [s.status for s in result.steps] == ["ok", status, "not_computed"]
    assert result.steps[1].message
    assert result.steps[1].properties is None


@pytest.mark.parametrize(
    "step", [NormalShockStep(), ExpansionStep(0.1), ConicalShockStep(0.1)]
)
def test_subsonic_rejection(step: FlowStep) -> None:
    assert (
        solve_flow_sequence(FlowState(10000, 300, 0.5), [step]).steps[0].status
        == "invalid_input"
    )


def test_total_out_of_range_does_not_discard_static() -> None:
    initial = FlowState(10000, 1500, 3)
    result = solve_flow_sequence(
        initial,
        [ExpansionStep(0.02), IsentropicStep(2), NormalShockStep()],
        AIR_HARMONIC_OSCILLATOR,
    )
    assert result.initial.total_temperature is None
    assert result.initial.total_pressure is None
    assert result.initial.warning
    assert result.steps[0].status == "ok"
    assert result.steps[0].total_pressure_recovery == 1
    assert result.steps[1].status == "out_of_range"
    assert result.steps[2].status == "not_computed"


def test_sonic_area_and_subsonic_branch() -> None:
    from aerophysics.isentropic import area_ratio

    initial = FlowState(10000, 300, 2)
    sonic = solve_flow_sequence(
        initial,
        [IsentropicStep(1 / float(area_ratio(2)), "area_ratio", MachBranch.SUBSONIC)],
    )
    assert sonic.complete
    assert exit_properties(sonic, 0).state.mach == pytest.approx(1)
    subsonic = solve_flow_sequence(
        initial, [IsentropicStep(1, "area_ratio", MachBranch.SUBSONIC)]
    )
    assert subsonic.complete
    assert exit_properties(subsonic, 0).state.mach < 1
    assert exit_properties(subsonic, 0).mass_flux == pytest.approx(
        subsonic.initial.mass_flux
    )


@pytest.mark.parametrize("error_type", [ExpansionConvergenceError, RuntimeError])
def test_convergence_failure(
    monkeypatch: pytest.MonkeyPatch, error_type: type[RuntimeError]
) -> None:
    def fail(*args: object, **kwargs: object) -> None:
        raise error_type("test failure")

    monkeypatch.setattr("aerophysics.flow_sequence.prandtl_meyer_expansion", fail)
    assert (
        solve_flow_sequence(FlowState(1, 300, 2), [ExpansionStep(0.1)]).steps[0].status
        == "convergence_error"
    )


def test_unsupported_step() -> None:
    with pytest.raises(TypeError):
        solve_flow_sequence(FlowState(1, 300, 2), [cast(FlowStep, object())])


def test_valid_shock_temperature_round_trip_at_upper_bound() -> None:
    initial = FlowState(25_000.0, 404.40881763527057, 4.850018713943831)
    result = solve_flow_sequence(initial, [NormalShockStep()], AIR_HARMONIC_OSCILLATOR)
    assert result.complete
    assert exit_properties(result, 0).state.temperature <= 2000.0
    assert exit_properties(result, 0).state.temperature == pytest.approx(2000, abs=2e-9)
    beyond = FlowState(initial.pressure, initial.temperature, initial.mach + 1e-9)
    invalid = solve_flow_sequence(
        beyond, [NormalShockStep(), ExpansionStep(0.1)], AIR_HARMONIC_OSCILLATOR
    )
    assert [step.status for step in invalid.steps] == ["out_of_range", "not_computed"]


@pytest.mark.parametrize("boundary,direction", [(400.0, -math.inf), (2000.0, math.inf)])
def test_restore_only_one_ulp_of_solver_output(
    boundary: float, direction: float
) -> None:
    from aerophysics._temperature import restore_static_temperature

    one = math.nextafter(boundary, direction)
    two = math.nextafter(one, direction)
    gas = AIR_HARMONIC_OSCILLATOR
    assert restore_static_temperature(one, gas) == boundary
    assert restore_static_temperature(two, gas) == two
    # Even one ULP outside is invalid for user inputs, not just two ULP.
    for temperature in (one, two):
        with pytest.raises(ModelRangeError):
            flow_properties(FlowState(25000, temperature, 1), gas)


@pytest.mark.parametrize("mach", [0.0, 2.0, 10.0])
def test_unbounded_oscillator_total_state_and_empty_sequence(mach: float) -> None:
    import warnings

    from aerophysics.gas import PerfectGas
    from aerophysics.real_gas import HarmonicOscillatorGas

    gas = HarmonicOscillatorGas(287.05, 1.4)
    initial = FlowState(100_000.0, 300.0, mach)
    expected = flow_properties(initial, PerfectGas(287.05, 1.4))
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        actual = flow_properties(initial, gas)
        sequence = solve_flow_sequence(initial, [], gas)
    assert actual.total_temperature == pytest.approx(300 * (1 + 0.2 * mach**2))
    assert actual.total_pressure == pytest.approx(expected.total_pressure, rel=2e-12)
    assert actual.total_enthalpy == pytest.approx(expected.total_enthalpy)
    assert actual.warning == ""
    assert sequence.complete
    assert sequence.initial == actual


def test_finite_thermal_ceiling_retains_static_properties() -> None:
    from aerophysics.real_gas import HarmonicOscillatorGas

    gas = HarmonicOscillatorGas(287.05, 1.4, applicable_temperature_range=(200, 500))
    result = flow_properties(FlowState(100000, 300, 2), gas)
    assert result.total_temperature is None
    assert result.total_pressure is None
    assert result.state.temperature == 300
    assert result.warning


def test_unbounded_vibrational_gas_conserves_stagnation_enthalpy() -> None:
    from aerophysics.real_gas import HarmonicOscillatorGas, VibrationalMode

    gas = HarmonicOscillatorGas(287.05, 1.4, modes=(VibrationalMode(1, 3055.56),))
    state = FlowState(100000, 500, 3)
    result = flow_properties(state, gas)
    assert result.total_temperature is not None
    assert 500 < result.total_temperature < 1400
    assert float(gas.standard_enthalpy(result.total_temperature)) == pytest.approx(
        float(gas.standard_enthalpy(500)) + result.velocity**2 / 2, rel=2e-12
    )
