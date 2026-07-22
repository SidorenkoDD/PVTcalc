"""Контракт физических диагностик результата без запуска EOS."""

from calc_core.Utils.ResultDiagnostics import diagnose_flash_result
from calc_core.VLE.FlashResult import FlashResult, PhaseState


def _result(**changes) -> FlashResult:
    values = {
        "pressure": 100.0,
        "temperature": 350.0,
        "vapor": PhaseState(
            0.4, {"A": 0.75, "B": 0.25},
            {"molar_volume": 10.0, "molar_density": 0.1, "density": 5.0},
        ),
        "liquid": PhaseState(
            0.6, {"A": 1 / 3, "B": 2 / 3},
            {"molar_volume": 2.0, "molar_density": 0.5, "density": 40.0},
        ),
        "is_two_phase": True,
    }
    values.update(changes)
    return FlashResult(**values)


def test_valid_flash_has_ok_quality() -> None:
    result = _result()
    feed = {"A": 0.5, "B": 0.5}

    diagnostics = diagnose_flash_result(result, feed)

    assert diagnostics.status == "ok"
    assert diagnostics.warnings == ()


def test_physical_problems_are_warnings_not_exceptions() -> None:
    result = _result(
        vapor=PhaseState(
            1.2, {"A": 0.8, "B": 0.4},
            {"molar_volume": -1.0, "density": float("nan")},
        ),
    )

    diagnostics = diagnose_flash_result(result, {"A": 0.5, "B": 0.5})
    codes = {item.code for item in diagnostics.warnings}

    assert diagnostics.status == "warning"
    assert {"phase_fraction_bounds", "phase_fraction_sum", "composition_sum",
            "non_positive_property", "non_finite", "material_balance"} <= codes


def test_diagnostics_round_trip_is_json_shaped() -> None:
    diagnostics = diagnose_flash_result(
        _result(pressure=float("inf")), {"A": 0.5, "B": 0.5},
    )

    restored = type(diagnostics).from_dict(diagnostics.to_dict())

    assert restored == diagnostics
    assert restored.status == "warning"

