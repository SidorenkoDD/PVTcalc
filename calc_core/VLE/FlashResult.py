from dataclasses import dataclass
from typing import Dict, Any

@dataclass(frozen=True)
class PhaseState:
    """Состояние одной фазы"""
    mole_fraction: float  # От 0.0 до 1.0 (доля этой фазы в системе)
    composition: Any      # Состав фазы (ваш объект Composition или dict/list)
    properties: Dict[str, Any] # Словарь свойств (плотность, Z-фактор и т.д.)

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

    @property
    def liquid_composition(self) -> Any:
        """Удобный геттер для DLE: всегда возвращает состав жидкости для следующего шага"""
        return self.liquid.composition
