"""Tests for NASA-polynomial thermochemistry and frozen dry air."""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics import (
    AIR,
    AIR_NASA7,
    AIR_NASA9,
    IdealGasSpecies,
    NASA7Polynomial,
    NASA9Polynomial,
    ThermallyPerfectGas,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.thermochemistry import (
    STANDARD_PRESSURE,
    UNIVERSAL_GAS_CONSTANT,
)


def test_nasa7_equations_and_array_shape() -> None:
    polynomial = NASA7Polynomial(
        (1.0, 10.0),
        ((1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),),
    )
    temperature = np.asarray([[2.0, 3.0]])

    cp = polynomial.cp_over_r(temperature)
    enthalpy = polynomial.h_over_rt(2.0)
    entropy = polynomial.s_over_r(2.0)

    assert isinstance(cp, np.ndarray)
    assert cp.shape == (1, 2)
    assert_allclose(cp, [[129.0, 547.0]])
    assert enthalpy == pytest.approx(34.0)
    assert entropy == pytest.approx(np.log(2.0) + 4.0 + 6.0 + 32.0 / 3.0 + 20.0 + 7.0)


def test_nasa9_equations() -> None:
    polynomial = NASA9Polynomial(
        (1.0, 10.0),
        ((1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0),),
    )
    temperature = 2.0

    assert polynomial.cp_over_r(temperature) == pytest.approx(
        1.0 / 4.0 + 1.0 + 3.0 + 8.0 + 20.0 + 48.0 + 112.0
    )
    assert polynomial.h_over_rt(temperature) == pytest.approx(
        -1.0 / 4.0 + np.log(2.0) + 3.0 + 4.0 + 20.0 / 3.0 + 12.0 + 112.0 / 5.0 + 4.0
    )
    assert polynomial.s_over_r(temperature) == pytest.approx(
        -1.0 / 8.0 - 1.0 + 3.0 * np.log(2.0) + 8.0 + 10.0 + 16.0 + 28.0 + 9.0
    )


@pytest.mark.parametrize("polynomial_type", [NASA7Polynomial, NASA9Polynomial])
def test_shared_boundary_uses_lower_region(
    polynomial_type: type[NASA7Polynomial] | type[NASA9Polynomial],
) -> None:
    coefficient_count = 7 if polynomial_type is NASA7Polynomial else 9
    low = (1.0,) + (0.0,) * (coefficient_count - 1)
    high = (2.0,) + (0.0,) * (coefficient_count - 1)
    if polynomial_type is NASA9Polynomial:
        low = (0.0, 0.0, 1.0) + (0.0,) * 6
        high = (0.0, 0.0, 2.0) + (0.0,) * 6
    polynomial = polynomial_type((200.0, 1000.0, 6000.0), (low, high))

    assert polynomial.cp_over_r(1000.0) == pytest.approx(1.0)
    assert polynomial.cp_over_r(np.nextafter(1000.0, np.inf)) == pytest.approx(2.0)
    assert polynomial.cp_over_r(6000.0) == pytest.approx(2.0)


@pytest.mark.parametrize(
    ("polynomial_type", "count"),
    [(NASA7Polynomial, 7), (NASA9Polynomial, 9)],
)
def test_polynomial_range_error_and_explicit_extrapolation(
    polynomial_type: type[NASA7Polynomial] | type[NASA9Polynomial],
    count: int,
) -> None:
    coefficients = (0.0, 0.0, 1.0) + (0.0,) * (count - 3)
    if polynomial_type is NASA7Polynomial:
        coefficients = (1.0,) + (0.0,) * 6
    polynomial = polynomial_type((200.0, 6000.0), (coefficients,))

    with pytest.raises(ModelRangeError, match="200 K and 6000 K"):
        polynomial.cp_over_r(6001.0)
    with pytest.warns(ApplicabilityWarning, match="extrapolated") as caught:
        result = polynomial.cp_over_r([100.0, 7000.0], allow_extrapolation=True)

    assert len(caught) == 1
    assert_allclose(result, [1.0, 1.0])


@pytest.mark.parametrize(
    ("ranges", "rows", "match"),
    [
        ((200.0,), (), "at least one"),
        ((200.0, 1000.0, 6000.0), ((1.0,) * 7,), "one more"),
        ((200.0, 200.0), ((1.0,) * 7,), "strictly increasing"),
        ((-1.0, 200.0), ((1.0,) * 7,), "positive"),
        ((200.0, 6000.0), ((1.0,) * 6,), "7 values"),
        ((200.0, np.inf), ((1.0,) * 7,), "finite"),
    ],
)
def test_polynomial_rejects_invalid_configuration(
    ranges: tuple[float, ...],
    rows: tuple[tuple[float, ...], ...],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        NASA7Polynomial(ranges, rows)


def test_polynomial_rejects_non_numeric_configuration() -> None:
    with pytest.raises(ValueError, match="real numeric"):
        NASA7Polynomial(
            (200.0, 6000.0),
            (("coefficient",) * 7,),  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan, np.inf, "hot"])
def test_polynomial_rejects_invalid_temperature(temperature: object) -> None:
    polynomial = NASA7Polynomial((200.0, 6000.0), ((1.0,) + (0.0,) * 6,))
    with pytest.raises(ValueError):
        polynomial.cp_over_r(temperature)  # type: ignore[arg-type]


def test_ideal_gas_species_standard_properties() -> None:
    polynomial = NASA7Polynomial(
        (200.0, 6000.0),
        ((3.5, 0.0, 0.0, 0.0, 0.0, 10.0, 2.0),),
    )
    species = IdealGasSpecies("X2", 0.028, polynomial)

    assert species.temperature_range == (200.0, 6000.0)
    assert species.specific_gas_constant == pytest.approx(
        UNIVERSAL_GAS_CONSTANT / 0.028
    )
    assert species.standard_molar_cp(300.0) == pytest.approx(
        3.5 * UNIVERSAL_GAS_CONSTANT
    )
    assert species.standard_molar_enthalpy(300.0) == pytest.approx(
        UNIVERSAL_GAS_CONSTANT * (3.5 * 300.0 + 10.0)
    )
    assert species.standard_molar_entropy(300.0) == pytest.approx(
        UNIVERSAL_GAS_CONSTANT * (3.5 * np.log(300.0) + 2.0)
    )


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"name": "", "molar_mass": 0.028}, ValueError),
        ({"name": "X", "molar_mass": 0.0}, ValueError),
        ({"name": "X", "molar_mass": np.inf}, ValueError),
        ({"name": "X", "molar_mass": 0.028, "reference_pressure": 0.0}, ValueError),
    ],
)
def test_ideal_gas_species_rejects_invalid_values(
    kwargs: dict[str, object], error: type[Exception]
) -> None:
    polynomial = NASA7Polynomial((200.0, 6000.0), ((1.0,) + (0.0,) * 6,))
    with pytest.raises(error):
        IdealGasSpecies(thermo=polynomial, **kwargs)  # type: ignore[arg-type]


def test_ideal_gas_species_rejects_invalid_polynomial_type() -> None:
    with pytest.raises(TypeError, match="NASA7Polynomial or NASA9Polynomial"):
        IdealGasSpecies("X", 0.028, object())  # type: ignore[arg-type]


def test_mixture_validates_composition_and_common_range() -> None:
    low = IdealGasSpecies(
        "low",
        0.02,
        NASA7Polynomial((200.0, 500.0), ((3.5,) + (0.0,) * 6,)),
    )
    high = IdealGasSpecies(
        "high",
        0.03,
        NASA7Polynomial((600.0, 1000.0), ((3.5,) + (0.0,) * 6,)),
    )

    with pytest.raises(ValueError, match="at least one"):
        ThermallyPerfectGas((), ())
    with pytest.raises(ValueError, match="same length"):
        ThermallyPerfectGas((low,), (0.5, 0.5))
    with pytest.raises(TypeError, match="IdealGasSpecies"):
        ThermallyPerfectGas((object(),), (1.0,))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="sum to one"):
        ThermallyPerfectGas((low,), (0.5,))
    with pytest.raises(ValueError, match="real numeric"):
        ThermallyPerfectGas((low,), ("all",))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="greater than zero"):
        ThermallyPerfectGas((low, high), (1.0, 0.0))
    with pytest.raises(ValueError, match="unique"):
        ThermallyPerfectGas((low, low), (0.5, 0.5))
    with pytest.raises(ValueError, match="common temperature"):
        ThermallyPerfectGas((low, high), (0.5, 0.5))


def test_built_in_air_composition_and_reference_values() -> None:
    assert AIR_NASA7.temperature_range == (200.0, 6000.0)
    assert AIR_NASA9.temperature_range == (200.0, 6000.0)
    assert AIR_NASA9.molar_mass == pytest.approx(0.028964766130783925)
    assert sum(AIR_NASA9.mass_fractions) == pytest.approx(1.0)
    assert AIR_NASA9.specific_gas_constant == pytest.approx(287.054367385917)
    assert AIR_NASA9.specific_gas_constant == pytest.approx(
        AIR.specific_gas_constant, rel=2e-5
    )

    temperatures = [200.0, 298.15, 1000.0, 3500.0, 6000.0]
    assert_allclose(
        AIR_NASA7.cp(temperatures),
        [
            1003.0871649585773,
            1004.7389659559424,
            1140.6749303994145,
            1308.4746668402383,
            1357.2479041100664,
        ],
        rtol=2e-13,
    )
    assert_allclose(
        AIR_NASA9.cp(temperatures),
        [
            1002.4061557658488,
            1004.7389671485392,
            1141.0326884939145,
            1309.3206938128276,
            1360.6655257996658,
        ],
        rtol=2e-13,
    )


@pytest.mark.parametrize("air", [AIR_NASA7, AIR_NASA9])
def test_mixture_thermodynamic_identities(air: ThermallyPerfectGas) -> None:
    temperature = np.asarray([[300.0, 1000.0], [2000.0, 6000.0]])
    cp = np.asarray(air.cp(temperature))
    cv = np.asarray(air.cv(temperature))
    gamma = np.asarray(air.heat_capacity_ratio(temperature))
    enthalpy = np.asarray(air.standard_enthalpy(temperature))
    internal_energy = np.asarray(air.standard_internal_energy(temperature))
    speed_of_sound = np.asarray(air.speed_of_sound(temperature))

    assert cp.shape == temperature.shape
    assert_allclose(cp - cv, air.specific_gas_constant)
    assert_allclose(gamma, cp / cv)
    assert_allclose(
        internal_energy,
        enthalpy - air.specific_gas_constant * temperature,
    )
    assert_allclose(
        speed_of_sound**2,
        gamma * air.specific_gas_constant * temperature,
    )


def test_sensible_properties_use_explicit_reference() -> None:
    temperature = [298.15, 1000.0]
    enthalpy = np.asarray(AIR_NASA9.sensible_enthalpy(temperature))
    internal_energy = np.asarray(AIR_NASA9.sensible_internal_energy(temperature))

    assert enthalpy[0] == pytest.approx(0.0, abs=1e-10)
    assert internal_energy[0] == pytest.approx(0.0, abs=1e-10)
    assert_allclose(
        enthalpy - internal_energy,
        AIR_NASA9.specific_gas_constant * (np.asarray(temperature) - temperature[0]),
    )
    assert AIR_NASA9.sensible_enthalpy(
        500.0, reference_temperature=500.0
    ) == pytest.approx(0.0, abs=1e-10)
    with pytest.raises(ValueError, match="finite scalar"):
        AIR_NASA9.sensible_enthalpy(500.0, reference_temperature=[300.0, 400.0])  # type: ignore[arg-type]


def test_entropy_includes_mixing_and_pressure() -> None:
    temperature = 1000.0
    entropy_at_reference = AIR_NASA9.entropy(temperature, STANDARD_PRESSURE)
    entropy_at_double_pressure = AIR_NASA9.entropy(temperature, 2.0 * STANDARD_PRESSURE)
    expected_difference = -AIR_NASA9.specific_gas_constant * np.log(2.0)

    assert entropy_at_double_pressure - entropy_at_reference == pytest.approx(
        expected_difference
    )

    standard_molar_entropy = sum(
        fraction * species.standard_molar_entropy(temperature)
        for species, fraction in zip(
            AIR_NASA9.species, AIR_NASA9.mole_fractions, strict=True
        )
    )
    mixing = -UNIVERSAL_GAS_CONSTANT * sum(
        fraction * np.log(fraction) for fraction in AIR_NASA9.mole_fractions
    )
    assert entropy_at_reference == pytest.approx(
        (standard_molar_entropy + mixing) / AIR_NASA9.molar_mass
    )


def test_entropy_broadcasting_and_invalid_pressure() -> None:
    entropy = AIR_NASA9.entropy(300.0, [[STANDARD_PRESSURE], [2e5]])
    assert isinstance(entropy, np.ndarray)
    assert entropy.shape == (2, 1)

    with pytest.raises(ValueError, match="greater than zero"):
        AIR_NASA9.entropy(300.0, 0.0)
    with pytest.raises(ValueError, match="broadcast-compatible"):
        AIR_NASA9.entropy([300.0, 400.0], [1e5, 2e5, 3e5])


def test_mixture_extrapolation_warns_once_per_public_call() -> None:
    with pytest.raises(ModelRangeError):
        AIR_NASA9.cp(7000.0)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        value = AIR_NASA9.heat_capacity_ratio([100.0, 7000.0], allow_extrapolation=True)

    assert np.all(np.isfinite(value))
    extrapolation_warnings = [
        warning
        for warning in caught
        if issubclass(warning.category, ApplicabilityWarning)
    ]
    assert len(extrapolation_warnings) == 1


@pytest.mark.parametrize(
    ("method", "args"),
    [
        ("sensible_enthalpy", (300.0, np.nan)),
        ("sensible_internal_energy", (300.0, -1.0)),
        ("standard_enthalpy", (np.nan,)),
        ("speed_of_sound", (0.0,)),
    ],
)
def test_mixture_rejects_invalid_temperature(
    method: str, args: tuple[object, ...]
) -> None:
    with pytest.raises(ValueError):
        getattr(AIR_NASA9, method)(*args)
