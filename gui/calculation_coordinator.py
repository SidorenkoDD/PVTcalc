"""Фреймворк-независимая оркестрация фоновых расчётов GUI."""

from __future__ import annotations

import logging
import threading
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

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
    done: threading.Event = field(default_factory=threading.Event)


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
              operation: Callable[[], Any]) -> CalculationJob:
        if self._active is not None:
            raise RuntimeError("Another calculation is already running")
        job = CalculationJob(node_ref=node_ref, label=label)
        self._active = job

        def worker() -> None:
            try:
                job.result = operation()
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
        return True

    def take_finished(self) -> CalculationJob | None:
        """Возвращает завершённую задачу и освобождает координатор."""
        job = self._active
        if job is None or not job.done.is_set():
            return None
        self._active = None
        return job
