"""Tests for GUI unit, configuration, and table helpers."""

import json

import numpy as np
import pytest

from aerophysics.detached_shock import DetachedShockGeometry
from aerophysics.gui.adapters import Row, detached_shock_shape
from aerophysics.gui.config import (
    CONFIG_SCHEMA_VERSION,
    ConfigurationError,
    dump_configuration,
    load_configuration,
    make_configuration,
    validate_configuration,
)
from aerophysics.gui.tables import (
    columns_for,
    detached_shock_shape_csv,
    display_rows,
    rows_to_csv,
)
from aerophysics.gui.units import (
    UnitPreferences,
    from_si,
    inverse_length_unit,
    selected_unit,
    to_si,
)


@pytest.mark.parametrize(
    ("kind", "unit", "values"),
    [
        ("length", "ft", [0.0, 3.0, 10.0]),
        ("length", "mm", [0.0, 3.0, 10.0]),
        ("length", "in", [0.0, 3.0, 10.0]),
        ("area", "ft²", [0.0, 3.0, 10.0]),
        ("area", "in²", [0.0, 3.0, 10.0]),
        ("speed", "kt", [0.0, 100.0]),
        ("speed", "ft/s", [0.0, 100.0]),
        ("pressure", "psi", [0.0, 14.7]),
        ("pressure", "psf", [0.0, 2116.0]),
        ("pressure", "kPa", [0.0, 101.325]),
        ("pressure", "hPa", [0.0, 1013.25]),
        ("density", "slug/ft³", [0.0, 0.00237]),
        ("density", "lbm/ft³", [0.0, 0.0765]),
        ("force", "lbf", [0.0, 100.0]),
        ("inverse_length", "1/mm", [0.0, 100.0]),
        ("inverse_length", "1/ft", [0.0, 100.0]),
        ("inverse_length", "1/in", [0.0, 100.0]),
        ("angle", "deg", [0.0, 10.0, 90.0]),
        ("temperature", "°C", [-40.0, 0.0, 100.0]),
        ("temperature", "°F", [-40.0, 32.0, 100.0]),
        ("temperature", "°R", [0.0, 491.67, 671.67]),
    ],
)
def test_display_unit_round_trip(kind: str, unit: str, values: list[float]) -> None:
    si = to_si(values, kind, unit)  # type: ignore[arg-type]
    assert np.asarray(from_si(si, kind, unit)) == pytest.approx(values)  # type: ignore[arg-type]


def test_scalar_units_and_preferences_validation() -> None:
    assert to_si(1.0, "length", "m") == 1.0
    preferences = UnitPreferences(
        length="ft",
        speed="kt",
        pressure="psi",
        temperature="°F",
        density="slug/ft³",
        angle="rad",
    )
    assert selected_unit("speed", preferences) == "kt"
    assert inverse_length_unit("ft") == "1/ft"
    assert UnitPreferences.from_dict(preferences.to_dict()) == preferences
    with pytest.raises(ValueError, match="unsupported length"):
        UnitPreferences(length="yard")
    with pytest.raises(ValueError, match="unsupported length"):
        inverse_length_unit("yard")
    with pytest.raises(ValueError, match="units must"):
        UnitPreferences.from_dict([])
    with pytest.raises(ValueError, match="missing"):
        UnitPreferences.from_dict({"length": "m"})
    with pytest.raises(ValueError, match="finite"):
        to_si(np.nan, "length", "m")
    with pytest.raises(ValueError, match="unsupported"):
        to_si(1.0, "pressure", "bar")


def test_configuration_round_trip() -> None:
    preferences = UnitPreferences(length="ft", angle="deg")
    configuration = make_configuration(
        calculator="flight",
        mode="sweep",
        inputs_si={
            "geometric_altitude": 1000.0,
            "motion": 0.8,
            "characteristic_length": None,
        },
        models={"motion_basis": "mach"},
        units=preferences,
        sweep_si={"field": "altitude", "start": 0.0, "stop": 2000.0, "points": 3},
    )
    serialized = dump_configuration(configuration)
    restored = load_configuration(serialized)
    assert restored == configuration
    assert restored["schema_version"] == CONFIG_SCHEMA_VERSION
    assert json.loads(serialized)["display_units"]["length"] == "ft"


def test_legacy_display_unit_configuration_uses_new_defaults() -> None:
    legacy = {
        "length": "ft",
        "speed": "kt",
        "pressure": "psi",
        "temperature": "°F",
        "density": "slug/ft³",
        "angle": "deg",
    }
    preferences = UnitPreferences.from_dict(legacy)
    assert preferences.length == "ft"
    assert preferences.area == "m²"
    assert preferences.force == "N"
    assert preferences.inverse_length == "1/m"


def test_detached_shock_configuration_table_and_shape_csv() -> None:
    configuration = make_configuration(
        calculator="detached_shock",
        mode="sweep",
        inputs_si={"upstream_mach": 4.0, "nose_radius": 0.1},
        models={
            "geometry": DetachedShockGeometry.AXISYMMETRIC_SPHERE.value,
            "model": "comparison",
        },
        units=UnitPreferences(length="ft"),
        sweep_si={"field": "upstream_mach", "start": 2.0, "stop": 8.0, "points": 5},
    )
    assert load_configuration(dump_configuration(configuration)) == configuration
    assert columns_for("detached_shock")[0].key == "upstream_mach"

    shape = detached_shock_shape(
        upstream_mach=4.0,
        nose_radius=0.1,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        points=3,
    )
    csv = detached_shock_shape_csv(shape, UnitPreferences(length="ft"))
    assert "x [ft],y [ft]" in csv
    assert len(csv.splitlines()) == 4


def test_configuration_round_trip_with_embedded_profile_arrays() -> None:
    configuration = make_configuration(
        calculator="protrusion_drag",
        mode="single",
        inputs_si={
            "drag_coefficient": 1.0,
            "height": 0.01,
            "base_width": 0.005,
            "edge_velocity": 100.0,
            "edge_density": 1.0,
            "boundary_layer_thickness": 0.02,
            "mach": None,
            "edge_temperature": None,
            "wall_temperature": None,
            "profile_height": [0.0, 0.01],
            "profile_velocity": [0.0, 100.0],
            "profile_density": [1.2, 1.0],
            "shape_height": None,
            "shape_width": None,
        },
        models={
            "profile_source": "csv",
            "shape": "rectangle",
            "compressible": False,
        },
        units=UnitPreferences(),
    )
    assert load_configuration(dump_configuration(configuration)) == configuration


@pytest.mark.parametrize(
    "value",
    [
        [],
        {},
        {
            "schema_version": 99,
            "calculator": "flight",
            "mode": "single",
            "inputs_si": {},
            "models": {},
            "display_units": UnitPreferences().to_dict(),
        },
    ],
)
def test_configuration_rejects_invalid_objects(value: object) -> None:
    with pytest.raises(ConfigurationError):
        validate_configuration(value)


def test_configuration_rejects_invalid_json_and_fields() -> None:
    with pytest.raises(ConfigurationError, match="valid JSON"):
        load_configuration("{")
    base = {
        "schema_version": 1,
        "calculator": "unknown",
        "mode": "single",
        "inputs_si": {},
        "models": {},
        "display_units": UnitPreferences().to_dict(),
    }
    with pytest.raises(ConfigurationError, match="calculator"):
        validate_configuration(base)
    base["calculator"] = "flight"
    base["mode"] = "other"
    with pytest.raises(ConfigurationError, match="mode"):
        validate_configuration(base)
    base["inputs_si"] = {
        "geometric_altitude": 0.0,
        "motion": 0.8,
        "characteristic_length": None,
    }
    base["models"] = {"motion_basis": "mach"}
    base["mode"] = "sweep"
    with pytest.raises(ConfigurationError, match="sweep_si"):
        validate_configuration(base)


def test_display_table_and_csv() -> None:
    rows = (
        {
            "upstream_mach": 2.0,
            "deflection_angle": np.deg2rad(10.0),
            "maximum_deflection_angle": np.deg2rad(23.0),
            "shock_angle": np.deg2rad(39.0),
            "downstream_mach": 1.64,
            "static_pressure_ratio": 1.7,
            "static_density_ratio": 1.45,
            "static_temperature_ratio": 1.17,
            "total_pressure_ratio": 0.98,
            "status": "ok",
            "message": "",
        },
    )
    table = display_rows("oblique_shock", rows, UnitPreferences())
    assert table[0]["偏向角 θ [deg]"] == pytest.approx(10.0)
    csv_text = rows_to_csv(table)
    assert csv_text.startswith("\ufeff")
    assert "偏向角 θ [deg]" in csv_text
    with pytest.raises(ValueError, match="empty"):
        rows_to_csv([])
    with pytest.raises(ValueError, match="unsupported"):
        columns_for("missing")


def test_display_table_converts_new_quantity_kinds() -> None:
    flight: tuple[Row, ...] = (
        {
            "reynolds_number_per_length": 1000.0,
        },
    )
    flight_table = display_rows(
        "flight", flight, UnitPreferences(inverse_length="1/ft")
    )
    assert flight_table[0]["Reynolds数/長さ [1/ft]"] == pytest.approx(304.8)

    isentropic: tuple[Row, ...] = ({"static_density": 1.225},)
    isentropic_table = display_rows(
        "isentropic", isentropic, UnitPreferences(density="lbm/ft³")
    )
    density_heading = next(
        column.heading(UnitPreferences(density="lbm/ft³"))
        for column in columns_for("isentropic")
        if column.key == "static_density"
    )
    assert isentropic_table[0][density_heading] == pytest.approx(0.07647, rel=1e-4)

    protrusion: tuple[Row, ...] = (
        {"direct_drag": 4.4482216152605, "frontal_area": 0.09290304},
    )
    protrusion_table = display_rows(
        "protrusion_drag",
        protrusion,
        UnitPreferences(force="lbf", area="ft²"),
    )
    assert protrusion_table[0]["直接抗力 D [lbf]"] == pytest.approx(1.0)
    assert protrusion_table[0]["前面面積 [ft²]"] == pytest.approx(1.0)


def test_conical_shock_table_converts_angles() -> None:
    rows = (
        {
            "upstream_mach": 2.0,
            "cone_half_angle": np.deg2rad(10.0),
            "maximum_cone_half_angle": np.deg2rad(40.0),
            "shock_angle": np.deg2rad(31.0),
            "post_shock_mach": 1.95,
            "surface_mach": 1.83,
            "surface_pressure_ratio": 1.29,
            "surface_density_ratio": 1.20,
            "surface_temperature_ratio": 1.08,
            "total_pressure_ratio": 0.999,
            "status": "ok",
            "message": "",
        },
    )
    table = display_rows("conical_shock", rows, UnitPreferences())
    assert table[0]["円錐半頂角 θc [deg]"] == pytest.approx(10.0)
    assert table[0]["衝撃波角 β [deg]"] == pytest.approx(31.0)


_VALID_SINGLE_PAYLOADS: dict[str, tuple[dict[str, object], dict[str, object]]] = {
    "flight": (
        {
            "geometric_altitude": 1000.0,
            "motion": 0.8,
            "characteristic_length": None,
        },
        {"motion_basis": "mach"},
    ),
    "oblique_shock": (
        {"upstream_mach": 2.0, "deflection_angle": 0.1},
        {"branch": "weak"},
    ),
    "conical_shock": (
        {"upstream_mach": 2.0, "cone_half_angle": 0.1},
        {},
    ),
    "boundary_layer": (
        {
            "distance": 1.0,
            "edge_velocity": 100.0,
            "edge_density": 1.0,
            "edge_dynamic_viscosity": 1.0e-5,
            "transition_reynolds": None,
            "mach": None,
            "edge_temperature": 288.15,
            "wall_temperature": None,
        },
        {
            "source": "manual",
            "regime": "turbulent",
            "turbulent_correlation": "schlichting",
            "compressibility_correction": "none",
        },
    ),
    "isentropic": (
        {
            "input_value": 2.0,
            "total_pressure": None,
            "total_temperature": 300.0,
        },
        {
            "input_basis": "mach",
            "branch": "subsonic",
            "gas_model": "AIR",
            "with_mass_flux": False,
            "allow_extrapolation": True,
        },
    ),
    "normal_shock": ({"upstream_mach": 2.0}, {}),
    "expansion": ({"upstream_mach": 2.0, "turn_angle": 0.1}, {}),
    "detached_shock": (
        {"upstream_mach": 4.0, "nose_radius": 0.1},
        {"geometry": "axisymmetric_sphere", "model": "comparison"},
    ),
    "boundary_layer_profile": (
        {
            "edge_velocity": 300.0,
            "edge_density": 1.0,
            "edge_temperature": 300.0,
            "boundary_layer_thickness": 0.05,
            "wall_shear_stress": 85.0,
            "wall_temperature": None,
            "wake_parameter": None,
            "points": 257,
        },
        {
            "source": "manual",
            "transformation": "compare",
            "temperature_velocity_relation": "generalized_reynolds_analogy",
        },
    ),
    "protrusion_drag": (
        {
            "drag_coefficient": 1.0,
            "height": 0.01,
            "base_width": 0.005,
            "edge_velocity": 100.0,
            "edge_density": 1.0,
            "boundary_layer_thickness": 0.02,
            "mach": None,
            "edge_temperature": None,
            "wall_temperature": None,
            "profile_height": None,
            "profile_velocity": None,
            "profile_density": None,
            "shape_height": None,
            "shape_width": None,
        },
        {
            "profile_source": "power_law",
            "shape": "rectangle",
            "compressible": False,
        },
    ),
    "thermochemistry": (
        {
            "temperature": 300.0,
            "pressure": 101_325.0,
            "reference_temperature": 298.15,
        },
        {"selection": "compare", "allow_extrapolation": False},
    ),
    "viscosity": (
        {"temperature": 1000.0},
        {"selection": "compare", "allow_extrapolation": False},
    ),
}


@pytest.mark.parametrize("calculator", sorted(_VALID_SINGLE_PAYLOADS))
def test_calculator_tables_and_configuration_schemas(calculator: str) -> None:
    assert columns_for(calculator)
    inputs_si, models = _VALID_SINGLE_PAYLOADS[calculator]
    configuration = make_configuration(
        calculator=calculator,
        mode="single",
        inputs_si=inputs_si,
        models=models,
        units=UnitPreferences(),
    )
    assert configuration["calculator"] == calculator


def test_viscosity_table_converts_temperature_and_preserves_fixed_units() -> None:
    rows: tuple[Row, ...] = (
        {
            "model": "Keyes",
            "temperature": 300.0,
            "dynamic_viscosity": 1.8519327e-5,
            "relative_difference": 0.320,
            "status": "ok",
            "message": "",
        },
    )
    table = display_rows("viscosity", rows, UnitPreferences(temperature="°F"))
    assert table[0]["温度 T [°F]"] == pytest.approx(80.33)
    assert table[0]["粘性係数 μ [Pa·s]"] == pytest.approx(1.8519327e-5)
    assert table[0]["Sutherland基準相対差 [%]"] == pytest.approx(0.320)
