"""Фреймворк-независимая оркестрация фоновых расчётов GUI."""

from __future__ import annotations

import inspect
import logging
import threading
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from calc_core.Utils.Cancellation import (
    CalculationCancelled,
    CancellationToken,
)
from gui.app_state import NodeRef

logger = logging.getLogger(__name__)


@dataclass
class CalculationJob:
    """Одна фоновая задача, адресованная полным ``NodeRef``."""

    node_ref: NodeRef
    label: str
    result: Any = None
    error: Optional[str] = None
    cancelled: bool = False
    progress: float = 0.0
    progress_message: str = ""
    cancellation_token: CancellationToken = field(default_factory=CancellationToken)
    done: threading.Event = field(default_factory=threading.Event)

    def report_progress(self, fraction: float, message: str = "") -> None:
        """Обновляет отображаемый прогресс фоновой задачи."""
        self.progress = max(0.0, min(1.0, float(fraction)))
        self.progress_message = message


class CalculationCoordinator:
    """Запускает не более одной задачи и не зависит от DearPyGui."""

    def __init__(self):
        self._active: CalculationJob | None = None

    @property
    def active(self) -> CalculationJob | None:
        return self._active

    @property
    def busy(self) -> bool:
        return self._active is not None

    def start(self, node_ref: NodeRef, label: str,
              operation: Callable[..., Any]) -> CalculationJob:
        if self._active is not None:
            raise RuntimeError("Another calculation is already running")
        job = CalculationJob(node_ref=node_ref, label=label)
        self._active = job

        def worker() -> None:
            try:
                # Старые zero-argument callbacks остаются совместимыми. Новые
                # операции принимают (CancellationToken, ProgressCallback).
                parameters = inspect.signature(operation).parameters.values()
                positional = [p for p in parameters if p.kind in (
                    inspect.Parameter.POSITIONAL_ONLY,
                    inspect.Parameter.POSITIONAL_OR_KEYWORD,
                )]
                accepts_varargs = any(
                    p.kind is inspect.Parameter.VAR_POSITIONAL
                    for p in parameters
                )
                if accepts_varargs or len(positional) >= 2:
                    job.result = operation(
                        job.cancellation_token, job.report_progress,
                    )
                elif len(positional) == 1:
                    job.result = operation(job.cancellation_token)
                else:
                    job.result = operation()
            except CalculationCancelled:
                job.cancelled = True
            except Exception as exc:  # noqa: BLE001 — передать ошибку в UI
                job.error = str(exc)
                logger.exception("Фоновый расчёт '%s' не удался", label)
            finally:
                job.done.set()

        threading.Thread(target=worker, daemon=True,
                         name=f"pvtcalc-{label}").start()
        return job

    def cancel(self) -> bool:
        """Запрашивает отмену; непрерываемый движок закончит работу в фоне."""
        if self._active is None:
            return False
        self._active.cancelled = True
        self._active.cancellation_token.cancel()
        return True

    def take_finished(self) -> CalculationJob | None:
        """Возвращает завершённую задачу и освобождает координатор."""
        job = self._active
        if job is None or not job.done.is_set():
            return None
        self._active = None
        return job
