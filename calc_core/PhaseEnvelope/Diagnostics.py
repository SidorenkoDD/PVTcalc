"""Общие статусы качества расчёта фазовой огибающей.

Численные методы по-прежнему возвращают привычную таблицу с ``NaN`` в
отсутствующих точках. Подробная причина каждого ``NaN`` прикладывается к
``DataFrame.attrs['diagnostics']`` и не меняет публичные колонки таблицы.
Это временный компактный контракт до унификации ``SaturationPoint`` в R1.3.
"""

from __future__ import annotations

from typing import TypedDict

import pandas as pd

OK = "ok"
NOT_PRESENT = "not_present"
NOT_SEARCHED = "not_searched"
SOLVER_NONCONVERGENCE = "solver_nonconvergence"
STABILITY_NONCONVERGENCE = "stability_nonconvergence"
VERIFICATION_FAILED = "verification_failed"
OUT_OF_RANGE = "out_of_range"
SKIPPED_AFTER_FAILURE = "skipped_after_failure"
SLOW_POINT_LIMIT = "slow_point_limit"

FAILURE_STATUSES = {
    SOLVER_NONCONVERGENCE,
    STABILITY_NONCONVERGENCE,
    VERIFICATION_FAILED,
    OUT_OF_RANGE,
    SKIPPED_AFTER_FAILURE,
    SLOW_POINT_LIMIT,
}


class EnvelopeDiagnostic(TypedDict):
    """Диагностика двух возможных границ фаз при одной температуре."""

    Temp_C: float
    primary_status: str
    secondary_status: str


def attach_diagnostics(
    frame: pd.DataFrame,
    diagnostics: list[EnvelopeDiagnostic],
    *,
    method: str,
    fallback: str | None,
) -> pd.DataFrame:
    """Прикладывает JSON-совместимую диагностику, не меняя колонки таблицы."""
    rows = [dict(row) for row in diagnostics]
    failed_temperatures = sorted({
        float(row["Temp_C"])
        for row in rows
        if row["primary_status"] in FAILURE_STATUSES
        or row["secondary_status"] in FAILURE_STATUSES
    })
    frame.attrs["diagnostics"] = {
        "method": method,
        "partial": bool(failed_temperatures),
        "failed_temperatures": failed_temperatures,
        "fallback": fallback,
        "rows": rows,
    }
    return frame
