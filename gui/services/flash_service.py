"""
Запуск одного флэш-расчёта (P, T) — фреймворк-независимо.

Тонкая адаптация `CompositionalModel.flash`: передаёт давление в барах и
температуру в °C, оставляя JSON-сериализацию слепка результата этому модулю.
Публичный фасад обрабатывает и двухфазный, и однофазный случай. Не импортирует
DearPyGui; оркестрация потока/прогресса/отмены — во View.

Перед созданием фасада сервис делает нормализованную глубокую копию состава,
чтобы округлённые сохранённые доли и общий объект редактора не попадали внутрь
расчёта.
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
from calc_core.Utils.ResultDiagnostics import ResultDiagnostics
from calc_core.VLE.FlashResult import FlashResult, PhaseState

logger = logging.getLogger(__name__)

# Строки таблицы результата: ключ в PhaseState.properties -> подпись.
# `mole_fraction` берётся из самого PhaseState, не из properties (см. View).
PHASE_PROPERTY_ROWS: list[tuple[str, str]] = [
    ("molecular_ weight", "Molar mass, g/mol"),  # ключ с пробелом — как в движке
    ("molar_volume", "Molar volume"),
    ("molar_density", "Molar density"),
    ("density", "Density"),
    ("z", "Z-factor"),
    ("viscosity", "Viscosity"),
]


def run_flash(composition: Composition, p_bar: float, t_celsius: float) -> FlashResult:
    """
    Считает флэш для состава при заданных давлении (бар) и температуре (°C).

    Returns
    -------
    FlashResult
        `pressure`, `temperature` (K), `vapor`/`liquid` (`PhaseState`),
        `is_two_phase`.
    """
    logger.info("Флэш: P=%s бар, T=%s °C", p_bar, t_celsius)
    # Сохранённые инженерные составы часто округлены до 5-6 знаков; фасад
    # принимает уже нормированный Composition, поэтому нормализуем копию здесь.
    work = composition.new_composition(composition.composition, deep_copy=True)
    return CompositionalModel(work).flash(p_bar, t_celsius)


def _as_float(value):
    """Приводит значение свойства к float (numpy/строки), иначе None — для JSON."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def snapshot_flash_result(result: FlashResult) -> dict:
    """
    Сериализуемый слепок результата флэша для сессии (только отображаемые
    числа: доли фаз + `properties`, без объекта состава). См. `restore`.
    """
    def phase(p: PhaseState) -> dict:
        comp = p.composition if isinstance(p.composition, dict) else None
        return {
            "mole_fraction": _as_float(p.mole_fraction),
            "properties": {k: _as_float(v) for k, v in p.properties.items()},
            "composition": ({k: _as_float(v) for k, v in comp.items()}
                            if comp else None),
        }

    return {
        "is_two_phase": bool(result.is_two_phase),
        "phase_type": result.phase_type,
        "diagnostics": result.diagnostics.to_dict(),
        "pressure": _as_float(result.pressure),
        "temperature": _as_float(result.temperature),
        "vapor": phase(result.vapor),
        "liquid": phase(result.liquid),
    }


def restore_flash_result(snap: dict) -> FlashResult:
    """Восстанавливает `FlashResult` из слепка (состав фаз не хранится → None)."""
    def phase(d: dict) -> PhaseState:
        return PhaseState(mole_fraction=d.get("mole_fraction"),
                          composition=d.get("composition"),
                          properties=d.get("properties", {}))

    return FlashResult(
        pressure=snap.get("pressure"),
        temperature=snap.get("temperature"),
        vapor=phase(snap.get("vapor", {})),
        liquid=phase(snap.get("liquid", {})),
        is_two_phase=bool(snap.get("is_two_phase")),
        phase_type=snap.get("phase_type"),
        diagnostics=ResultDiagnostics.from_dict(snap.get("diagnostics")),
    )
