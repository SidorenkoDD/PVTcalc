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
from calc_core.VLE.FlashResult import FlashResult

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
