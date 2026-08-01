"""Tests for gas transport-property models."""

import math
from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.gas import AIR_CONDUCTIVITY as LEGACY_AIR_CONDUCTIVITY
from aerophysics.gas import AIR_VISCOSITY as LEGACY_AIR_VISCOSITY
from aerophysics.gas import SutherlandModel as LegacySutherlandModel
from aerophysics.gas import USSAConductivityModel as LegacyUSSAConductivityModel
from aerophysics.transport import (
    AIR_BLOTTNER_VISCOSITY,
    AIR_CONDUCTIVITY,
    AIR_KEYES_VISCOSITY,
    AIR_VISCOSITY,
    BlottnerModel,
    DynamicViscosityModel,
    KeyesModel,
    SutherlandModel,
    USSAConductivityModel,
    WilkeMixtureViscosityModel,
)


def _blottner_value(temperature: float, a: float, b: float, c: float) -> float:
    logarithm = math.log(temperature)
    return 0.1 * math.exp((a * logarithm + b) * logarithm + c)


def _wilke_value(
    viscosities: tuple[float, ...],
    molar_masses: tuple[float, ...],
    mole_fractions: tuple[float, ...],
) -> float:
    result = 0.0
    for i, viscosity_i in enumerate(viscosities):
        denominator = 0.0
        for j, viscosity_j in enumerate(viscosities):
            phi = (
                1.0
                + math.sqrt(viscosity_i / viscosity_j)
                * (molar_masses[j] / molar_masses[i]) ** 0.25
            ) ** 2 / math.sqrt(8.0 * (1.0 + molar_masses[i] / molar_masses[j]))
            denominator += mole_fractions[j] * phi
        result += mole_fractions[i] * viscosity_i / denominator
    return result


def test_legacy_gas_transport_exports_are_identity_aliases() -> None:
    assert LEGACY_AIR_VISCOSITY is AIR_VISCOSITY
    assert LEGACY_AIR_CONDUCTIVITY is AIR_CONDUCTIVITY
    assert LegacySutherlandModel is SutherlandModel
    assert LegacyUSSAConductivityModel is USSAConductivityModel


def test_keyes_reference_value_and_array_shape() -> None:
    temperature = 300.0
    expected = (
        1.488e-6
        * temperature**1.5
        / (temperature + 122.1 * 10.0 ** (-5.0 / temperature))
    )
    assert AIR_KEYES_VISCOSITY.dynamic_viscosity(temperature) == pytest.approx(
        expected
    )
    result = AIR_KEYES_VISCOSITY.dynamic_viscosity([[100.0, 300.0], [500.0, 1000.0]])
    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 2)
    assert result.dtype == np.float64


def test_blottner_reference_value_and_array_shape() -> None:
    model = BlottnerModel(0.0268142, 0.3177838, -11.3155513)
    assert model.dynamic_viscosity(1000.0) == pytest.approx(
        _blottner_value(1000.0, model.a, model.b, model.c)
    )
    result = model.dynamic_viscosity([[1000.0, 2000.0], [5000.0, 30_000.0]])
    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 2)
    assert result.dtype == np.float64


@pytest.mark.parametrize(
    ("model", "lower", "upper"),
    [
        (AIR_KEYES_VISCOSITY, 79.0, 1845.0),
        (BlottnerModel(0.0268142, 0.3177838, -11.3155513), 1000.0, 30_000.0),
    ],
)
def test_correlation_boundaries_and_extrapolation_warning(
    model: KeyesModel | BlottnerModel, lower: float, upper: float
) -> None:
    model.dynamic_viscosity([lower, upper])
    with pytest.warns(ApplicabilityWarning, match="correlation was extrapolated"):
        result = model.dynamic_viscosity([0.5 * lower, 2.0 * upper])
    assert np.all(np.isfinite(result))


def test_wilke_single_component_and_identical_components() -> None:
    single = WilkeMixtureViscosityModel((AIR_VISCOSITY,), (0.028,), (1.0,))
    expected = AIR_VISCOSITY.dynamic_viscosity([300.0, 1000.0])
    assert_allclose(single.dynamic_viscosity([300.0, 1000.0]), expected)

    identical = WilkeMixtureViscosityModel(
        (AIR_VISCOSITY, AIR_VISCOSITY),
        (0.028, 0.028),
        (0.3, 0.7),
    )
    assert_allclose(identical.dynamic_viscosity([300.0, 1000.0]), expected)


def test_wilke_permutation_invariance_and_dry_air_reference() -> None:
    models = AIR_BLOTTNER_VISCOSITY.component_models
    masses = AIR_BLOTTNER_VISCOSITY.molar_masses
    fractions = AIR_BLOTTNER_VISCOSITY.mole_fractions
    viscosities = tuple(
        float(model.dynamic_viscosity(1000.0)) for model in models
    )
    expected = _wilke_value(viscosities, masses, fractions)
    assert AIR_BLOTTNER_VISCOSITY.dynamic_viscosity(1000.0) == pytest.approx(
        expected
    )
    assert expected == pytest.approx(4.137574698616173e-5)

    reversed_model = WilkeMixtureViscosityModel(
        tuple(reversed(models)),
        tuple(reversed(masses)),
        tuple(reversed(fractions)),
    )
    assert reversed_model.dynamic_viscosity(1000.0) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ((0.0, 122.1, 5.0), "coefficient"),
        ((1.0, -1.0, 5.0), "additive_temperature"),
        ((1.0, 122.1, np.inf), "exponential_temperature"),
    ],
)
def test_keyes_rejects_invalid_coefficients(
    arguments: tuple[float, float, float], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        KeyesModel(*arguments)


@pytest.mark.parametrize(
    "arguments",
    [
        (np.nan, 0.0, 0.0),
        (0.0, np.inf, 0.0),
        (0.0, 0.0, -np.inf),
    ],
)
def test_blottner_rejects_nonfinite_coefficients(
    arguments: tuple[float, float, float],
) -> None:
    with pytest.raises(ValueError, match="must be finite"):
        BlottnerModel(*arguments)


@pytest.mark.parametrize(
    "temperature_range",
    [(0.0, 1000.0), (1000.0, 1000.0), (1000.0, np.inf)],
)
def test_models_reject_invalid_temperature_range(
    temperature_range: tuple[float, float],
) -> None:
    with pytest.raises(ValueError, match="temperature_range"):
        BlottnerModel(0.0, 0.0, 0.0, temperature_range)


def test_model_rejects_malformed_temperature_range() -> None:
    malformed = cast(tuple[float, float], (1000.0,))
    with pytest.raises(ValueError, match="two real numeric values"):
        BlottnerModel(0.0, 0.0, 0.0, malformed)


@pytest.mark.parametrize(
    ("models", "masses", "fractions", "message"),
    [
        ((), (), (), "at least one"),
        ((AIR_VISCOSITY,), (0.028, 0.032), (1.0,), "same length"),
        ((AIR_VISCOSITY,), (0.0,), (1.0,), "molar_masses"),
        ((AIR_VISCOSITY,), (0.028,), (0.0,), "mole_fractions"),
        (
            (AIR_VISCOSITY, AIR_VISCOSITY),
            (0.028, 0.032),
            (0.4, 0.5),
            "sum to one",
        ),
    ],
)
def test_wilke_rejects_invalid_definition(
    models: tuple[SutherlandModel, ...],
    masses: tuple[float, ...],
    fractions: tuple[float, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        WilkeMixtureViscosityModel(models, masses, fractions)


def test_wilke_rejects_invalid_component_and_non_numeric_data() -> None:
    invalid_model = cast(DynamicViscosityModel, object())
    with pytest.raises(TypeError, match="dynamic-viscosity models"):
        WilkeMixtureViscosityModel((invalid_model,), (0.028,), (1.0,))

    invalid_masses = cast(tuple[float, ...], ("heavy",))
    with pytest.raises(ValueError, match="real numeric values"):
        WilkeMixtureViscosityModel(
            (AIR_VISCOSITY,), invalid_masses, (1.0,)
        )


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan, np.inf])
def test_new_viscosity_models_reject_invalid_temperature(temperature: float) -> None:
    with pytest.raises(ValueError):
        AIR_KEYES_VISCOSITY.dynamic_viscosity(temperature)
    with pytest.raises(ValueError):
        AIR_BLOTTNER_VISCOSITY.dynamic_viscosity(temperature)
