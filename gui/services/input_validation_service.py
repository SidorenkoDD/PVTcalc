"""Fast, GUI-only validation for calculation input forms.

The functions in this module mirror only the inexpensive engineering checks
needed to prevent an invalid job from being submitted.  Calculation services
keep their own authoritative validation, so this layer never changes solver
behaviour or the calculation core.
"""

from __future__ import annotations

import math
from typing import Any

_ABSOLUTE_ZERO_C = -273.15


def _number(value: object, label: str) -> tuple[float | None, str | None]:
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        return None, f"{label} must be a finite number"
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None, f"{label} must be a finite number"
    if not math.isfinite(result):
        return None, f"{label} must be a finite number"
    return result, None


def parse_number_list(value: object, label: str) -> tuple[list[float], str | None]:
    """Parses comma-, whitespace- or line-separated finite numbers."""
    if not isinstance(value, str):
        return [], f"{label} are required"
    parts = value.replace("\n", " ").replace(",", " ").split()
    if not parts:
        return [], f"{label} are required"
    numbers: list[float] = []
    for part in parts:
        number, error = _number(part, label)
        if error is not None or number is None:
            return [], f"{label} must contain only finite numbers"
        numbers.append(number)
    return numbers, None


def validate_flash_inputs(pressure: object, temperature_c: object) -> dict[str, str]:
    """Checks the pressure/temperature pair entered for one flash run."""
    errors: dict[str, str] = {}
    p, error = _number(pressure, "Pressure")
    if error is not None:
        errors["P"] = error
    elif p is not None and p <= 0:
        errors["P"] = "Pressure must be greater than 0 bar"
    t, error = _number(temperature_c, "Temperature")
    if error is not None:
        errors["T"] = error
    elif t is not None and t <= _ABSOLUTE_ZERO_C:
        errors["T"] = "Temperature must be above absolute zero"
    return errors


def validate_experiment_inputs(kind: str, pressure_text: object, temperature_c: object,
                               p_res: object | None = None,
                               stage_temperature_text: object | None = None,
                               ) -> tuple[dict[str, Any], dict[str, str]]:
    """Validates and returns normalized CCE/DLE/Separator parameters."""
    errors: dict[str, str] = {}
    pressures, error = parse_number_list(pressure_text, "Pressure stages")
    if error is not None:
        errors["pressures"] = error
    elif len(pressures) < 2:
        errors["pressures"] = "At least 2 pressure stages are required"
    elif any(value <= 0 for value in pressures):
        errors["pressures"] = "Pressure stages must be greater than 0 bar"
    elif any(first <= second for first, second in zip(pressures, pressures[1:])):
        errors["pressures"] = "Pressure stages must be strictly descending"

    temperature, error = _number(temperature_c, "Temperature")
    if error is not None:
        errors["T_c"] = error
    elif temperature is not None and temperature <= _ABSOLUTE_ZERO_C:
        errors["T_c"] = "Temperature must be above absolute zero"

    params: dict[str, Any] = {"pressures": pressures, "T_c": temperature}
    if kind in {"dle", "separator"}:
        reservoir_pressure, error = _number(p_res, "P res")
        if error is not None:
            errors["P_res"] = error
        elif reservoir_pressure is not None and reservoir_pressure <= 0:
            errors["P_res"] = "P res must be greater than 0 bar"
        params["P_res"] = reservoir_pressure
    if kind == "separator":
        stage_temperatures, error = parse_number_list(
            stage_temperature_text, "Stage temperatures")
        if error is not None:
            errors["stage_temps"] = error
        elif len(stage_temperatures) != len(pressures):
            errors["stage_temps"] = "Stage temperature count must match pressure stages"
        elif any(value <= _ABSOLUTE_ZERO_C for value in stage_temperatures):
            errors["stage_temps"] = "Stage temperatures must be above absolute zero"
        params["stage_temps_c"] = stage_temperatures
    return params, errors


def validate_envelope_inputs(method: str, values: dict[str, object]) -> dict[str, str]:
    """Checks the active phase-envelope parameter group before submission."""
    errors: dict[str, str] = {}
    if method == "grid":
        for key, label in (("grid_t_max_c", "T max"),
                           ("grid_p_max_bar", "P max")):
            value, error = _number(values.get(key), label)
            if error is not None:
                errors[key] = error
            elif value is not None and value <= 0:
                errors[key] = f"{label} must be greater than 0"
        for key, label in (("grid_t_points", "T points"),
                           ("grid_p_points", "P points")):
            value, error = _number(values.get(key), label)
            if error is not None:
                errors[key] = error
            elif value is not None and value < 2:
                errors[key] = f"{label} must be at least 2"
        return errors

    t_min, min_error = _number(values.get("t_min_c"), "T min")
    t_max, max_error = _number(values.get("t_max_c"), "T max")
    step, step_error = _number(values.get("t_step_c"), "T step")
    pressure, pressure_error = _number(values.get("p_max_bar"), "P max")
    if min_error is not None:
        errors["t_min_c"] = min_error
    elif t_min is not None and t_min <= _ABSOLUTE_ZERO_C:
        errors["t_min_c"] = "T min must be above absolute zero"
    if max_error is not None:
        errors["t_max_c"] = max_error
    elif t_max is not None and t_max <= _ABSOLUTE_ZERO_C:
        errors["t_max_c"] = "T max must be above absolute zero"
    if t_min is not None and t_max is not None and t_max <= t_min:
        errors["t_max_c"] = "T max must be greater than T min"
    if step_error is not None:
        errors["t_step_c"] = step_error
    elif step is not None and step <= 0:
        errors["t_step_c"] = "T step must be greater than 0"
    if pressure_error is not None:
        errors["p_max_bar"] = pressure_error
    elif pressure is not None and pressure <= 0:
        errors["p_max_bar"] = "P max must be greater than 0 bar"
    return errors
