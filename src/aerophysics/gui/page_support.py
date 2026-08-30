"""Typed support for the Streamlit calculator page pipeline.

Imported configuration is normalized before pages see it.  These helpers keep
the remaining page boundary deterministic: absent or non-numeric defaults use
the caller's explicit fallback, while malformed stored results and result-row
metrics are rejected instead of being rendered as partially missing output.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, NamedTuple

import numpy as np
import streamlit as st

from aerophysics.gui.units import QuantityKind, from_si, to_si


class PageDataError(ValueError):
    """A stored page result does not satisfy the rendering contract."""


class ConfigurationDefaults(NamedTuple):
    """Normalized configuration sections consumed by calculator pages."""

    inputs: dict[str, Any]
    models: dict[str, Any]
    sweep: dict[str, Any]


def configuration_defaults(
    configuration: Mapping[str, object] | None,
) -> ConfigurationDefaults:
    """Return independent dictionaries for the three optional page sections.

    Missing or non-mapping sections are treated as empty.  Individual widgets
    must still provide an explicit fallback through :func:`numeric_default`.
    """

    if configuration is None:
        return ConfigurationDefaults({}, {}, {})

    def section(name: str) -> dict[str, Any]:
        value = configuration.get(name)
        return dict(value) if isinstance(value, Mapping) else {}

    return ConfigurationDefaults(
        section("inputs_si"),
        section("models"),
        section("sweep_si"),
    )


def display_value(value: float, kind: QuantityKind, unit: str) -> float:
    """Convert one SI page value to its selected display unit."""

    return float(from_si(value, kind, unit))


def si_value(value: float, kind: QuantityKind, unit: str) -> float:
    """Convert one page input from its selected display unit to SI."""

    return float(to_si(value, kind, unit))


def numeric_default(values: Mapping[str, object], key: str, default: float) -> float:
    """Read a numeric page default, falling back for missing or invalid data."""

    value = values.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return default
    return float(value)


def array_default(values: Mapping[str, object], key: str) -> np.ndarray | None:
    """Read a stored one-dimensional array, or return ``None`` when absent."""

    value = values.get(key)
    if not isinstance(value, list):
        return None
    try:
        result = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise PageDataError(f"{key} must be a one-dimensional numeric array") from error
    if result.ndim != 1:
        raise PageDataError(f"{key} must be a one-dimensional numeric array")
    return result


def session_payload[ResultT](
    key: str, result_type: type[ResultT]
) -> tuple[ResultT, dict[str, object]] | None:
    """Return a typed calculate-to-render payload from Streamlit state.

    Missing payloads are normal before the first calculation.  A value under
    the expected session key is not normal unless it is the exact two-part
    ``(result, configuration)`` contract, so malformed values raise
    :class:`PageDataError` rather than silently hiding stale output.
    """

    if key not in st.session_state:
        return None
    value = st.session_state[key]
    if not isinstance(value, tuple) or len(value) != 2:
        raise PageDataError(f"{key} must contain a result/configuration pair")
    result, configuration = value
    if not isinstance(result, result_type) or not isinstance(configuration, dict):
        raise PageDataError(f"{key} contains an invalid result/configuration pair")
    return result, configuration


def metric_value(row: Mapping[str, object], contains: str) -> float | None:
    """Return an available metric by exact heading or unambiguous fragment."""

    headings = (
        [contains]
        if contains in row
        else [heading for heading in row if contains in heading]
    )
    if len(headings) != 1:
        raise PageDataError(
            f"expected exactly one result field containing {contains!r}; "
            f"found {len(headings)}"
        )
    value = row[headings[0]]
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PageDataError(f"result field {headings[0]!r} must be numeric")
    return float(value)


def render_metric(row: Mapping[str, object], label: str, contains: str) -> None:
    """Render a validated numeric result-row metric."""

    value = metric_value(row, contains)
    st.metric(label, f"{value:.5g}" if value is not None else "—")


__all__ = [
    "ConfigurationDefaults",
    "PageDataError",
    "array_default",
    "configuration_defaults",
    "display_value",
    "metric_value",
    "numeric_default",
    "render_metric",
    "session_payload",
    "si_value",
]
