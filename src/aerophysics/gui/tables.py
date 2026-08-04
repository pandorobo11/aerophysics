"""Localized result tables and CSV export."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass

from aerophysics.gui.adapters import CellValue, Row
from aerophysics.gui.units import (
    QuantityKind,
    UnitPreferences,
    from_si,
    selected_unit,
)


@dataclass(frozen=True, slots=True)
class Column:
    key: str
    label: str
    kind: QuantityKind | None = None
    fixed_unit: str | None = None

    def heading(self, preferences: UnitPreferences) -> str:
        unit = (
            selected_unit(self.kind, preferences)
            if self.kind is not None
            else self.fixed_unit
        )
        return f"{self.label} [{unit}]" if unit else self.label

    def convert(self, value: CellValue, preferences: UnitPreferences) -> CellValue:
        if not isinstance(value, float) or self.kind is None:
            return value
        return float(from_si(value, self.kind, selected_unit(self.kind, preferences)))


FLIGHT_COLUMNS = (
    Column("geometric_altitude", "幾何高度 h", "length"),
    Column("temperature", "静温 T", "temperature"),
    Column("pressure", "静圧 p", "pressure"),
    Column("density", "密度 ρ", "density"),
    Column("speed_of_sound", "音速 a", "speed"),
    Column("mach", "Mach M"),
    Column("velocity", "速度 V", "speed"),
    Column("dynamic_pressure", "動圧 q", "pressure"),
    Column("reynolds_number_per_length", "Reynolds数/長さ", fixed_unit="1/m"),
    Column("reynolds_number", "Reynolds数 Re"),
    Column("total_temperature", "全温 T₀", "temperature"),
    Column("total_pressure", "全圧 p₀", "pressure"),
    Column("status", "status"),
    Column("message", "message"),
)

SHOCK_COLUMNS = (
    Column("upstream_mach", "上流 Mach M₁"),
    Column("deflection_angle", "偏向角 θ", "angle"),
    Column("maximum_deflection_angle", "最大付着偏向角 θmax", "angle"),
    Column("shock_angle", "衝撃波角 β", "angle"),
    Column("downstream_mach", "下流 Mach M₂"),
    Column("static_pressure_ratio", "p₂/p₁"),
    Column("static_density_ratio", "ρ₂/ρ₁"),
    Column("static_temperature_ratio", "T₂/T₁"),
    Column("total_pressure_ratio", "p₀₂/p₀₁"),
    Column("status", "status"),
    Column("message", "message"),
)

CONICAL_SHOCK_COLUMNS = (
    Column("upstream_mach", "上流 Mach M∞"),
    Column("cone_half_angle", "円錐半頂角 θc", "angle"),
    Column("maximum_cone_half_angle", "最大付着半頂角 θc,max", "angle"),
    Column("shock_angle", "衝撃波角 β", "angle"),
    Column("post_shock_mach", "衝撃波直後 Mach M₂"),
    Column("surface_mach", "表面 Mach Mₛ"),
    Column("surface_pressure_ratio", "pₛ/p∞"),
    Column("surface_density_ratio", "ρₛ/ρ∞"),
    Column("surface_temperature_ratio", "Tₛ/T∞"),
    Column("total_pressure_ratio", "p₀₂/p₀∞"),
    Column("status", "status"),
    Column("message", "message"),
)

ISENTROPIC_COLUMNS = (
    Column("input_value", "入力値"),
    Column("input_basis", "入力基準"),
    Column("mach", "Mach M"),
    Column("total_temperature_ratio", "T₀/T"),
    Column("total_pressure_ratio", "p₀/p"),
    Column("total_density_ratio", "ρ₀/ρ"),
    Column("area_ratio", "A/A*"),
    Column("mass_flow_parameter", "質量流量パラメータ"),
    Column("mass_flux", "質量流束", fixed_unit="kg/(m²·s)"),
    Column("choked_mass_flux", "チョーク質量流束", fixed_unit="kg/(m²·s)"),
    Column("critical_temperature_ratio", "T₀/T*"),
    Column("critical_pressure_ratio", "p₀/p*"),
    Column("critical_density_ratio", "ρ₀/ρ*"),
    Column("status", "status"),
    Column("message", "message"),
)

NORMAL_SHOCK_COLUMNS = (
    Column("upstream_mach", "上流 Mach M₁"),
    Column("downstream_mach", "下流 Mach M₂"),
    Column("static_pressure_ratio", "p₂/p₁"),
    Column("static_density_ratio", "ρ₂/ρ₁"),
    Column("static_temperature_ratio", "T₂/T₁"),
    Column("total_pressure_ratio", "p₀₂/p₀₁"),
    Column("pitot_pressure_ratio", "p₀₂/p₁"),
    Column("status", "status"),
    Column("message", "message"),
)

EXPANSION_COLUMNS = (
    Column("upstream_mach", "上流 Mach M₁"),
    Column("turn_angle", "膨張角 θ", "angle"),
    Column("maximum_turn_angle", "最大膨張角", "angle"),
    Column("downstream_mach", "下流 Mach M₂"),
    Column("upstream_prandtl_meyer_angle", "ν₁", "angle"),
    Column("downstream_prandtl_meyer_angle", "ν₂", "angle"),
    Column("static_pressure_ratio", "p₂/p₁"),
    Column("static_density_ratio", "ρ₂/ρ₁"),
    Column("static_temperature_ratio", "T₂/T₁"),
    Column("status", "status"),
    Column("message", "message"),
)

BOUNDARY_LAYER_COLUMNS = (
    Column("distance", "前縁からの距離 x", "length"),
    Column("reynolds_number", "Reynolds数 Re_x"),
    Column("effective_reynolds_number", "有効 Reynolds数 Re_eff"),
    Column("boundary_layer_thickness", "境界層厚さ δ₉₉", "length"),
    Column("displacement_thickness", "排除厚さ δ*", "length"),
    Column("momentum_thickness", "運動量厚さ θ", "length"),
    Column("local_skin_friction_coefficient", "局所摩擦係数 C_f"),
    Column("average_skin_friction_coefficient", "平均摩擦係数 C̄_f"),
    Column("wall_shear_stress", "壁面せん断応力 τ_w", "pressure"),
    Column("drag_per_unit_width", "単位幅抗力 D′", fixed_unit="N/m"),
    Column("recovery_temperature", "回復温度 T_r", "temperature"),
    Column("wall_temperature", "壁温 T_w", "temperature"),
    Column("status", "status"),
    Column("message", "message"),
)

BOUNDARY_LAYER_PROFILE_COLUMNS = (
    Column("model", "変換モデル"),
    Column("wall_distance", "壁面距離 y", "length"),
    Column("outer_coordinate", "y/δ₉₉"),
    Column("wall_distance_plus", "y⁺"),
    Column("transformed_wall_coordinate", "変換壁座標"),
    Column("velocity", "速度 U", "speed"),
    Column("velocity_ratio", "U/U_e"),
    Column("velocity_plus", "U⁺"),
    Column("transformed_velocity_plus", "変換速度 U⁺"),
    Column("temperature", "温度 T", "temperature"),
    Column("density", "密度 ρ", "density"),
    Column("dynamic_viscosity", "粘性係数 μ", fixed_unit="Pa·s"),
    Column("local_mach_number", "局所 Mach M"),
    Column("dynamic_pressure", "動圧 q", "pressure"),
    Column("friction_velocity", "摩擦速度 u_τ", "speed"),
    Column("friction_reynolds_number", "摩擦 Reynolds数 Re_τ"),
    Column("recovery_temperature", "回復温度 T_r", "temperature"),
    Column("wall_temperature", "壁温 T_w", "temperature"),
    Column("wake_parameter", "wake parameter Π"),
    Column("displacement_thickness", "排除厚さ δ*", "length"),
    Column("momentum_thickness", "運動量厚さ θ", "length"),
    Column("shape_factor", "形状係数 H"),
    Column("local_skin_friction_coefficient", "局所摩擦係数 C_f"),
    Column("status", "status"),
    Column("message", "message"),
)

PROTRUSION_COLUMNS = (
    Column("drag_coefficient", "抗力係数 C_D"),
    Column("height", "突起高さ h", "length"),
    Column("base_width", "代表幅 b₀", "length"),
    Column("boundary_layer_thickness", "境界層厚さ δ", "length"),
    Column("shape", "投影形状"),
    Column("mach", "Mach M_e"),
    Column("direct_drag", "直接抗力 D", fixed_unit="N"),
    Column("effective_dynamic_pressure", "実効動圧 q_eff", "pressure"),
    Column("shielding_factor", "遮蔽係数"),
    Column("frontal_area", "前面面積", fixed_unit="m²"),
    Column("edge_dynamic_pressure", "外縁動圧 q_e", "pressure"),
    Column("height_to_boundary_layer_thickness", "h/δ"),
    Column("profile", "プロファイル"),
    Column("compressibility_applied", "圧縮性適用"),
    Column("status", "status"),
    Column("message", "message"),
)

THERMOCHEMISTRY_COLUMNS = (
    Column("model", "モデル"),
    Column("temperature", "温度 T", "temperature"),
    Column("pressure", "圧力 p", "pressure"),
    Column("reference_temperature", "基準温度 T_ref", "temperature"),
    Column("molar_mass", "モル質量", fixed_unit="kg/mol"),
    Column("specific_gas_constant", "比気体定数 R", fixed_unit="J/(kg·K)"),
    Column("cp", "定圧比熱 c_p", fixed_unit="J/(kg·K)"),
    Column("cv", "定容比熱 c_v", fixed_unit="J/(kg·K)"),
    Column("heat_capacity_ratio", "比熱比 γ"),
    Column("speed_of_sound", "音速 a", "speed"),
    Column("standard_enthalpy", "標準エンタルピー h°", fixed_unit="J/kg"),
    Column("standard_internal_energy", "標準内部エネルギー u°", fixed_unit="J/kg"),
    Column("sensible_enthalpy", "顕熱エンタルピー Δh", fixed_unit="J/kg"),
    Column("sensible_internal_energy", "顕熱内部エネルギー Δu", fixed_unit="J/kg"),
    Column("entropy", "エントロピー s", fixed_unit="J/(kg·K)"),
    Column("status", "status"),
    Column("message", "message"),
)


def columns_for(calculator: str) -> tuple[Column, ...]:
    """Return ordered table columns for a calculator."""
    choices = {
        "conical_shock": CONICAL_SHOCK_COLUMNS,
        "flight": FLIGHT_COLUMNS,
        "isentropic": ISENTROPIC_COLUMNS,
        "normal_shock": NORMAL_SHOCK_COLUMNS,
        "oblique_shock": SHOCK_COLUMNS,
        "expansion": EXPANSION_COLUMNS,
        "boundary_layer": BOUNDARY_LAYER_COLUMNS,
        "boundary_layer_profile": BOUNDARY_LAYER_PROFILE_COLUMNS,
        "protrusion_drag": PROTRUSION_COLUMNS,
        "thermochemistry": THERMOCHEMISTRY_COLUMNS,
    }
    try:
        return choices[calculator]
    except KeyError as error:
        raise ValueError(f"unsupported calculator: {calculator}") from error


def display_rows(
    calculator: str, rows: tuple[Row, ...], preferences: UnitPreferences
) -> list[dict[str, CellValue]]:
    """Convert SI rows into localized, unit-labelled display rows."""
    columns = columns_for(calculator)
    return [
        {
            column.heading(preferences): column.convert(
                row.get(column.key), preferences
            )
            for column in columns
        }
        for row in rows
    ]


def rows_to_csv(rows: list[dict[str, CellValue]]) -> str:
    """Encode display rows as Excel-friendly UTF-8 CSV."""
    if not rows:
        raise ValueError("cannot export an empty table")
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return "\ufeff" + output.getvalue()
