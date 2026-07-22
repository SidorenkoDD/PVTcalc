"""Cooperative cancellation and progress reporting for long calculations."""

from __future__ import annotations

import threading
from collections.abc import Callable

from calc_core.Utils.Errors import PVTCalcError

ProgressCallback = Callable[[float, str], None]


class CalculationCancelled(PVTCalcError):
    """Расчёт остановлен в безопасной точке по запросу вызывающей стороны."""


class CancellationToken:
    """Потокобезопасный флаг кооперативной отмены."""

    def __init__(self) -> None:
        self._event = threading.Event()

    def cancel(self) -> None:
        self._event.set()

    @property
    def is_cancelled(self) -> bool:
        return self._event.is_set()

    def throw_if_cancelled(self) -> None:
        if self.is_cancelled:
            raise CalculationCancelled("Расчёт отменён пользователем.")


def report_progress(
    callback: ProgressCallback | None,
    fraction: float,
    message: str = "",
) -> None:
    """Безопасно передаёт progress callback, если вызывающий его указал."""
    if callback is None:
        return
    callback(max(0.0, min(1.0, float(fraction))), message)
