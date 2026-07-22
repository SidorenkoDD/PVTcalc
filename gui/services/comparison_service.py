"""Подготовка данных для сравнения однотипных PVT-экспериментов.

Модуль не зависит от DearPyGui: проверяет совместимость результатов, строит
общие серии и отклонения относительно первого выбранного расчёта.
"""

from __future__ import annotations

import math
from typing import Any

from gui.services.experiment_service import series_for_plot


class IncompatibleComparisonError(ValueError):
    """Выбранные узлы нельзя содержательно сравнить."""


def _finite(value: object) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)))


def _grid(result: dict) -> list[float]:
    columns = result.get("columns", [])
    x_name = result.get("x")
    if x_name not in columns:
        return []
    x_index = columns.index(x_name)
    return [float(row[x_index]) for row in result.get("rows", [])
            if isinstance(row, list) and x_index < len(row) and _finite(row[x_index])]


def _same_grid(left: list[float], right: list[float], tolerance: float) -> bool:
    return len(left) == len(right) and all(
        math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)
        for a, b in zip(left, right)
    )


def _value_map(result: dict, column: str) -> dict[float, float]:
    xs, ys = series_for_plot(result, column)
    return {round(float(x), 9): float(y) for x, y in zip(xs, ys)}


def _same_parameter(left: object, right: object, tolerance: float) -> bool:
    if (isinstance(left, (int, float)) and not isinstance(left, bool)
            and isinstance(right, (int, float)) and not isinstance(right, bool)
            and _finite(left) and _finite(right)):
        return math.isclose(float(left), float(right), rel_tol=tolerance,
                            abs_tol=tolerance)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(
            _same_parameter(a, b, tolerance) for a, b in zip(left, right)
        )
    return left == right


def build_experiment_comparison(
    members: list[dict[str, Any]],
    *,
    grid_tolerance: float = 1e-8,
) -> dict[str, Any]:
    """Возвращает JSON-совместимые серии и отклонения для ≥2 экспериментов.

    Каждый участник содержит ``id``, ``label``, ``kind``, ``result`` и
    необязательный ``lab_data``. Разные сетки допустимы: кривые строятся на
    собственных X, а отклонения считаются только в общих точках давления.
    """
    if len(members) < 2:
        raise IncompatibleComparisonError("Для сравнения нужны минимум два результата.")
    kinds = {member.get("kind") for member in members}
    if len(kinds) != 1 or None in kinds:
        raise IncompatibleComparisonError(
            "Сравнивать можно только эксперименты одного типа.",
        )
    if any(not isinstance(member.get("result"), dict) for member in members):
        raise IncompatibleComparisonError("У каждого участника должен быть результат.")

    results = [member["result"] for member in members]
    common_columns = list(results[0].get("plot_all", []))
    for result in results[1:]:
        available = set(result.get("plot_all", []))
        common_columns = [column for column in common_columns if column in available]
    if not common_columns:
        raise IncompatibleComparisonError("У результатов нет общих числовых кривых.")

    grids = [_grid(result) for result in results]
    grid_compatible = all(
        _same_grid(grids[0], grid, grid_tolerance) for grid in grids[1:]
    )
    warnings = [] if grid_compatible else [
        "Pressure grids differ; curves use their own grids and deviations "
        "are shown only at matching pressure points.",
    ]
    reference_params = members[0].get("params")
    reference_params = reference_params if isinstance(reference_params, dict) else {}
    differing_conditions = []
    for member in members[1:]:
        params = member.get("params")
        params = params if isinstance(params, dict) else {}
        different = [field for field in ("T_c", "P_res", "stage_temps_c")
                     if not _same_parameter(reference_params.get(field),
                                            params.get(field), grid_tolerance)]
        if different:
            differing_conditions.append(
                f"{member.get('label', member.get('id', 'run'))}: "
                + ", ".join(different),
            )
    if differing_conditions:
        warnings.append(
            "Experiment conditions differ from the reference ("
            + "; ".join(differing_conditions)
            + "). Treat curves as different experiments.",
        )
    stale = [str(member.get("label", member.get("id", "run")))
             for member in members if member.get("stale") is True]
    if stale:
        warnings.append("Stale result(s): " + ", ".join(stale) + ". Recalculate first.")

    series: dict[str, list[dict[str, Any]]] = {}
    deviations: dict[str, list[dict[str, Any]]] = {}
    for column in common_columns:
        column_series = []
        for member in members:
            xs, ys = series_for_plot(member["result"], column)
            lab = member.get("lab_data")
            lab_x, lab_y = series_for_plot(lab, column) if isinstance(lab, dict) else ([], [])
            column_series.append({
                "id": str(member.get("id", "")),
                "label": str(member.get("label", member.get("id", "run"))),
                "x": xs,
                "y": ys,
                "lab_x": lab_x,
                "lab_y": lab_y,
            })
        series[column] = column_series

        reference = _value_map(results[0], column)
        rows = []
        for member, result in zip(members[1:], results[1:]):
            candidate = _value_map(result, column)
            for pressure in sorted(reference.keys() & candidate.keys(), reverse=True):
                ref_value = reference[pressure]
                candidate_value = candidate[pressure]
                delta = candidate_value - ref_value
                relative = (100.0 * delta / ref_value) if ref_value != 0.0 else None
                rows.append({
                    "pressure": pressure,
                    "member": str(member.get("label", member.get("id", "run"))),
                    "reference": ref_value,
                    "candidate": candidate_value,
                    "absolute": abs(delta),
                    "relative_percent": relative,
                })
        deviations[column] = rows

    return {
        "kind": next(iter(kinds)),
        "reference": str(members[0].get("label", members[0].get("id", "run"))),
        "grid_compatible": grid_compatible,
        "warnings": warnings,
        "columns": common_columns,
        "series": series,
        "deviations": deviations,
    }
