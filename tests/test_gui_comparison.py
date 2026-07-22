"""Чистый сервис R3.2: общие кривые, сетки и отклонения."""

import pytest

from gui.services.comparison_service import (
    IncompatibleComparisonError,
    build_experiment_comparison,
)


def _result(pressures, values, *, kind="dle"):
    return {
        "kind": kind,
        "x": "pressure",
        "columns": ["pressure", "Bo", "Rs"],
        "rows": [[pressure, value, value * 100]
                 for pressure, value in zip(pressures, values)],
        "plot_all": ["Bo", "Rs"],
    }


def _member(member_id, pressures, values, *, kind="dle", lab_data=None,
            params=None, stale=False):
    return {
        "id": member_id,
        "label": member_id.upper(),
        "kind": kind,
        "result": _result(pressures, values, kind=kind),
        "lab_data": lab_data,
        "params": params or {"T_c": 100.0, "P_res": 400.0},
        "stale": stale,
    }


def test_same_grid_builds_curves_and_deviations() -> None:
    comparison = build_experiment_comparison([
        _member("base", [300, 200, 100], [1.3, 1.2, 1.1]),
        _member("tuned", [300, 200, 100], [1.4, 1.1, 1.1]),
    ])

    assert comparison["grid_compatible"] is True
    assert comparison["warnings"] == []
    assert comparison["columns"] == ["Bo", "Rs"]
    assert len(comparison["series"]["Bo"]) == 2
    first = comparison["deviations"]["Bo"][0]
    assert first["pressure"] == 300
    assert first["absolute"] == pytest.approx(0.1)
    assert first["relative_percent"] == pytest.approx(100 * 0.1 / 1.3)


def test_different_grids_warn_and_use_only_matching_points() -> None:
    comparison = build_experiment_comparison([
        _member("base", [300, 200, 100], [1.3, 1.2, 1.1]),
        _member("other", [300, 150, 100], [1.4, 1.15, 1.0]),
    ])

    assert comparison["grid_compatible"] is False
    assert comparison["warnings"]
    assert [row["pressure"] for row in comparison["deviations"]["Bo"]] == [300, 100]


def test_lab_measurements_are_returned_for_overlay() -> None:
    lab = {"x": "pressure", "columns": ["pressure", "Bo"],
           "rows": [[250, 1.25], [150, 1.15]]}
    comparison = build_experiment_comparison([
        _member("base", [300, 200], [1.3, 1.2], lab_data=lab),
        _member("other", [300, 200], [1.4, 1.1]),
    ])

    base_series = comparison["series"]["Bo"][0]
    assert base_series["lab_x"] == [150, 250]
    assert base_series["lab_y"] == [1.15, 1.25]


def test_different_experiment_types_are_rejected() -> None:
    with pytest.raises(IncompatibleComparisonError, match="одного типа"):
        build_experiment_comparison([
            _member("dle", [300, 200], [1.3, 1.2]),
            _member("cce", [300, 200], [1.0, 1.1], kind="cce"),
        ])


def test_different_conditions_and_stale_results_are_explicit() -> None:
    comparison = build_experiment_comparison([
        _member("base", [300, 200], [1.3, 1.2]),
        _member("tuned", [300, 200], [1.4, 1.1],
                params={"T_c": 110.0, "P_res": 400.0}, stale=True),
    ])

    assert any("conditions differ" in warning for warning in comparison["warnings"])
    assert any("Stale result" in warning for warning in comparison["warnings"])
