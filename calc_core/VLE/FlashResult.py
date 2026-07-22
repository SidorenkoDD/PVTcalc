from dataclasses import dataclass
from typing import Any

from calc_core.Utils.ResultDiagnostics import ResultDiagnostics


@dataclass(frozen=True)
class PhaseState:
    """Состояние одной фазы"""
    mole_fraction: float  # От 0.0 до 1.0 (доля этой фазы в системе)
    composition: Any      # Состав фазы (ваш объект Composition или dict/list)
    properties: dict[str, Any]  # Словарь свойств (плотность, Z-фактор и т.д.)


@dataclass(frozen=True)
class FlashResult:
    """
    Унифицированный результат флеша.
    ВСЕГДА содержит и vapor, и liquid.
    """
    pressure : float
    temperature : float
    vapor : PhaseState
    liquid : PhaseState
    is_two_phase : bool
    phase_type: str | None = None
    diagnostics: ResultDiagnostics = ResultDiagnostics()

    @property
    def quality_status(self) -> str:
        """``ok`` либо ``warning`` для отображения во внешнем интерфейсе."""
        return self.diagnostics.status

    @property
    def liquid_composition(self) -> Any:
        """Удобный геттер для DLE: всегда возвращает состав жидкости для следующего шага"""
        return self.liquid.composition
