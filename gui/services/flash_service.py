"""
Запуск одного флэш-расчёта (P, T) — фреймворк-независимо.

Тонкая обёртка над `calc_core.VLE.Flash` + `Conditions`: обрабатывает и
двухфазный, и однофазный случай (в отличие от `CompositionalModel.flash`,
который ждёт двухфазную точку). Ввод давления — в барах, температуры — в
°C (внутри `Conditions` конвертируется в Кельвины). Не импортирует
DearPyGui; оркестрация потока/прогресса/отмены — во View.

Побочный эффект (из движка): `Flash` перезаписывает `composition.T`
температурой точки и пересчитывает EOS-зависимые свойства состава.
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.Utils.Conditions import Conditions
from calc_core.VLE.Flash import Flash
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
    conditions = Conditions(p_bar, t_celsius)
    logger.info("Флэш: P=%s бар, T=%s °C", p_bar, t_celsius)
    return Flash(composition, conditions).calculate()


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
        return {
            "mole_fraction": _as_float(p.mole_fraction),
            "properties": {k: _as_float(v) for k, v in p.properties.items()},
        }

    return {
        "is_two_phase": bool(result.is_two_phase),
        "pressure": _as_float(result.pressure),
        "temperature": _as_float(result.temperature),
        "vapor": phase(result.vapor),
        "liquid": phase(result.liquid),
    }


def restore_flash_result(snap: dict) -> FlashResult:
    """Восстанавливает `FlashResult` из слепка (состав фаз не хранится → None)."""
    def phase(d: dict) -> PhaseState:
        return PhaseState(mole_fraction=d.get("mole_fraction"),
                          composition=None, properties=d.get("properties", {}))

    return FlashResult(
        pressure=snap.get("pressure"),
        temperature=snap.get("temperature"),
        vapor=phase(snap.get("vapor", {})),
        liquid=phase(snap.get("liquid", {})),
        is_two_phase=bool(snap.get("is_two_phase")),
    )
