"""Analytical relations and regression checks for heating correlations."""

import re
import warnings
from dataclasses import FrozenInstanceError, fields
from typing import Any, cast

import numpy as np
import pytest

from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.gas import AIR
from aerophysics.heating import (
    StagnationHeatingModel as Model,
)
from aerophysics.heating import (
    StagnationHeatingResult,
)
from aerophysics.heating import (
    sphere_stagnation_heat_flux as sphere,
)
from aerophysics.heating import (
    swept_cylinder_heat_flux as cylinder,
)

ARGS = (10.0, 226.5, 1197.0, 0.1)


def sphere_a(**kwargs: Any) -> StagnationHeatingResult:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        return sphere(*ARGS, **kwargs)


def test_reference_values() -> None:
    fay = sphere_a(wall_temperature=300)
    tauber = sphere_a(wall_temperature=300, model=Model.TAUBER)
    assert fay.freestream_density == pytest.approx(0.018410422066669437, rel=1e-12)
    assert fay.freestream_velocity == pytest.approx(3017.0271650438317, rel=1e-12)
    assert fay.stagnation_temperature == pytest.approx(4756.5, rel=1e-12)
    assert fay.stagnation_pressure == pytest.approx(154672.7111953016, rel=1e-10)
    assert fay.velocity_gradient == pytest.approx(16460.871554347686, rel=1e-10)
    assert fay.heat_flux == pytest.approx(1986435.0556857858, rel=1e-10)
    assert tauber.heat_flux == pytest.approx(2020351.9615749859, rel=1e-10)
    assert tauber.velocity_gradient is None
    assert 0.9 < cast(float, fay.heat_flux) / cast(float, tauber.heat_flux) < 1.1
    for angle, expected in [
        (0, 1424182.5302905634),
        (np.pi / 6, 1174325.8608920162),
        (np.pi / 3, 609841.1935253066),
    ]:
        result = cylinder(*ARGS, angle, wall_temperature=300)
        assert result.heat_flux == pytest.approx(expected, rel=1e-10)
    assert cylinder(
        *ARGS, np.pi / 6, wall_temperature=300
    ).recovery_enthalpy == pytest.approx(4573982.589592741, rel=1e-10)
    assert (
        abs(cast(float, cylinder(*ARGS, np.pi / 2, wall_temperature=300).heat_flux))
        < 1e-6
    )
    assert cylinder(*ARGS, 0, wall_temperature=300).heat_flux / cast(
        float, tauber.heat_flux
    ) == pytest.approx(1.29 / 1.83, rel=1e-12)


@pytest.mark.parametrize("model", list(Model))
def test_sphere_scaling_and_sign(model: Model) -> None:
    base = sphere_a(wall_temperature=300, model=model)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        large = sphere(*ARGS[:3], 0.4, wall_temperature=300, model=model)
    assert large.heat_flux == pytest.approx(base.heat_flux / 2, rel=1e-12)
    assert sphere_a(
        wall_temperature=base.stagnation_temperature, model=model
    ).heat_flux == pytest.approx(0, abs=1e-6)
    assert (
        sphere_a(
            wall_temperature=base.stagnation_temperature + 100, model=model
        ).heat_flux
        < 0
    )


def test_cylinder_scaling_and_sign() -> None:
    base = cylinder(*ARGS, np.pi / 6, wall_temperature=300)
    assert cylinder(
        *ARGS[:3], 0.4, np.pi / 6, wall_temperature=300
    ).heat_flux == pytest.approx(base.heat_flux / 2, rel=1e-12)
    wall = base.recovery_enthalpy / AIR.cp
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        assert cylinder(
            *ARGS, np.pi / 6, wall_temperature=wall
        ).heat_flux == pytest.approx(0, abs=1e-6)
        assert cylinder(*ARGS, np.pi / 6, wall_temperature=wall + 100).heat_flux < 0


@pytest.mark.parametrize("index,value", [(0, 1), (0, 0), (1, 0), (2, -1), (3, 0)])
@pytest.mark.parametrize("is_cylinder", [False, True])
def test_invalid_inputs(index: int, value: float, is_cylinder: bool) -> None:
    args: Any = list(ARGS)
    args[index] = value
    with pytest.raises(ValueError):
        if is_cylinder:
            cylinder(args[0], args[1], args[2], args[3], 0, wall_temperature=300)
        else:
            sphere(*args, wall_temperature=300)


@pytest.mark.parametrize("is_cylinder", [False, True])
def test_other_invalid_inputs(is_cylinder: bool) -> None:
    def call(*args: Any, **kw: Any) -> Any:
        if is_cylinder:
            return cylinder(args[0], args[1], args[2], args[3], 0, **kw)
        return sphere(*args, **kw)

    for wall in (0, float("nan")):
        with pytest.raises(ValueError):
            call(*ARGS, wall_temperature=wall)
    with pytest.raises(ValueError, match="heating inputs must be broadcastable"):
        call([2, 3], [200, 300, 400], *ARGS[2:], wall_temperature=300)


def test_model_prandtl_and_sweep_validation() -> None:
    with pytest.raises(ValueError):
        sphere(*ARGS, wall_temperature=300, model=cast(Any, "tauber"))
    for pr in (0, -1, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            sphere(*ARGS, wall_temperature=300, prandtl_number=pr)
    assert (
        sphere(
            *ARGS,
            wall_temperature=300,
            model=Model.TAUBER,
            prandtl_number=float("nan"),
            viscosity_model=cast(Any, None),
        ).heat_flux
        > 0
    )
    for angle in (-0.01, np.pi / 2 + 0.01):
        with pytest.raises(ValueError):
            cylinder(*ARGS, angle, wall_temperature=300)


def test_warnings() -> None:
    with pytest.warns(ApplicabilityWarning, match="Fay-Riddell") as caught:
        sphere(*ARGS, wall_temperature=300)
    assert len(caught) == 1
    with warnings.catch_warnings():
        warnings.simplefilter("error", ApplicabilityWarning)
        sphere(3, 220, 1197, 0.1, wall_temperature=300)
        sphere(*ARGS, wall_temperature=300, model=Model.TAUBER)
        cylinder(*ARGS, 0, wall_temperature=300)
    for call in (
        lambda: sphere(50, 300, 1197, 0.1, wall_temperature=300, model=Model.TAUBER),
        lambda: cylinder(50, 300, 1197, 0.1, 0, wall_temperature=300),
    ):
        with pytest.warns(ApplicabilityWarning, match="14 km/s") as caught:
            call()
        assert len(caught) == 1
    for call in (
        lambda: sphere(*ARGS, wall_temperature=4100, model=Model.TAUBER),
        lambda: cylinder(*ARGS, 0, wall_temperature=4100),
    ):
        with pytest.warns(ApplicabilityWarning, match="4000 K") as caught:
            call()
        assert len(caught) == 1
    for call in (
        lambda: sphere(50, 300, 1197, 0.1, wall_temperature=4100, model=Model.TAUBER),
        lambda: cylinder(50, 300, 1197, 0.1, 0, wall_temperature=4100),
    ):
        with pytest.warns(ApplicabilityWarning) as caught:
            call()
        assert len(caught) == 2


@pytest.mark.parametrize("is_cylinder", [False, True])
def test_tauber_high_speed_assumption(is_cylinder: bool) -> None:
    def evaluate(mach: Any) -> Any:
        if is_cylinder:
            return cylinder(mach, 226.5, 1197, 0.1, np.pi / 3, wall_temperature=300)
        return sphere(mach, 226.5, 1197, 0.1, wall_temperature=300, model=Model.TAUBER)

    for mach in (7.0, [7.0, 7.2]):
        with pytest.warns(
            ApplicabilityWarning, match=re.escape("h_inf << V^2/2")
        ) as caught:
            evaluate(mach)
        assert len(caught) == 1
    with warnings.catch_warnings():
        warnings.simplefilter("error", ApplicabilityWarning)
        evaluate(7.2)


@pytest.mark.parametrize("model", list(Model))
def test_arrays(model: Model) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        result = sphere(
            np.array([[3], [4]]),
            220,
            1197,
            [0.1, 0.4],
            wall_temperature=300,
            model=model,
        )
        cyl = cylinder(
            np.array([[3], [4]]), 220, 1197, [0.1, 0.4], np.pi / 6, wall_temperature=300
        )
        assert isinstance(result.heat_flux, np.ndarray)
        assert isinstance(cyl.heat_flux, np.ndarray)
        for i, mach in enumerate((3, 4)):
            for j, radius in enumerate((0.1, 0.4)):
                scalar_sphere = sphere(
                    mach, 220, 1197, radius, wall_temperature=300, model=model
                )
                scalar_cylinder = cylinder(
                    mach, 220, 1197, radius, np.pi / 6, wall_temperature=300
                )
                for array_result, scalar_result in (
                    (result, scalar_sphere),
                    (cyl, scalar_cylinder),
                ):
                    for field in fields(array_result):
                        array_value = getattr(array_result, field.name)
                        scalar_value = getattr(scalar_result, field.name)
                        if isinstance(array_value, np.ndarray):
                            assert isinstance(scalar_value, float)
                            assert array_value[i, j] == pytest.approx(scalar_value)
                        else:
                            assert array_value == scalar_value
    for obj in (result, cyl):
        for field in fields(obj):
            val = getattr(obj, field.name)
            if isinstance(val, np.ndarray):
                assert (
                    val.shape == (2, 2)
                    and val.dtype == np.float64
                    and not val.flags.writeable
                )
        with pytest.raises(FrozenInstanceError):
            obj.heat_flux = 0  # type: ignore[misc]
    assert isinstance(sphere_a(wall_temperature=300).heat_flux, float)
    assert isinstance(cylinder(*ARGS, 0, wall_temperature=300).heat_flux, float)
