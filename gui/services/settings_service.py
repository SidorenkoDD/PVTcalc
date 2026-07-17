"""
Настройки движка для панели Settings (пока «показать/хранить»).

Собирает актуальные значения констант, стандартных условий и критериев
сходимости из движка (`calc_core.Utils.Constants`, `Conditions.StandardConditions`)
и хранит пользовательские правки в `gui_settings.json`.

**Важно:** на этом этапе правки НЕ влияют на расчёт — движок продолжает
использовать свои значения из модулей (многие импортируются по значению при
загрузке модуля, поэтому надёжная проводка — отдельная задача). Панель
показывает текущие значения и запоминает изменённые на будущее.
"""

import json
import logging
from pathlib import Path

from calc_core.Utils import Constants
from calc_core.Utils.Conditions import StandardConditions

logger = logging.getLogger(__name__)

DEFAULT_SETTINGS_PATH = "gui_settings.json"

# Ключи стандартных условий (не константы модуля — читаются из StandardConditions).
_STD_P = "STD_P"
_STD_T = "STD_T"

# Схема панели: группы -> поля (key, подпись, формат отображения).
SCHEMA: list[dict] = [
    {
        "group": "Constants",
        "fields": [
            {"key": "CONSTANT_R", "label": "Universal gas constant R", "fmt": "%.6g"},
        ],
    },
    {
        "group": "Standard conditions",
        "fields": [
            {"key": _STD_P, "label": "Standard pressure, bar", "fmt": "%.6g"},
            {"key": _STD_T, "label": "Standard temperature, K", "fmt": "%.6g"},
        ],
    },
    {
        "group": "Convergence criteria (tolerances)",
        "fields": [
            {"key": "TOL_TWO_PHASE_STABILITY_CONVERGENCE",
             "label": "Stability test convergence", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION",
             "label": "Stability trivial solution", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE",
             "label": "Flash bisection (Rachford-Rice)", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_FLASH_CONVERGENCE",
             "label": "Flash fugacity convergence", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_FLASH_TRIVIAL_SOLUTION",
             "label": "Flash trivial solution", "fmt": "%.3g"},
            {"key": "TOL_SAT_PRESSURE",
             "label": "Saturation pressure convergence", "fmt": "%.3g"},
        ],
    },
]

# все ключи схемы по порядку (для сериализации/сравнения)
ALL_KEYS = [f["key"] for grp in SCHEMA for f in grp["fields"]]


def engine_defaults() -> dict:
    """Актуальные значения из движка (то, чем он реально считает сейчас)."""
    std = StandardConditions()
    values = {name: getattr(Constants, name)
              for name in ALL_KEYS if hasattr(Constants, name)}
    values[_STD_P] = std.p
    values[_STD_T] = std.t
    return values


def load_settings(path: str = DEFAULT_SETTINGS_PATH) -> dict:
    """
    Значения для панели: дефолты движка, поверх которых наложены сохранённые
    пользовательские правки (если файл есть и валиден).
    """
    values = engine_defaults()
    p = Path(path)
    if p.exists():
        try:
            saved = json.loads(p.read_text(encoding="utf-8"))
            for k in ALL_KEYS:
                if k in saved:
                    values[k] = float(saved[k])
        except (json.JSONDecodeError, TypeError, ValueError) as exc:
            logger.warning("Настройки %s повреждены (%s), беру дефолты движка", p, exc)
    return values


def save_settings(values: dict, path: str = DEFAULT_SETTINGS_PATH) -> None:
    """Пишет только известные ключи схемы (как float) в `gui_settings.json`."""
    out = {k: float(values[k]) for k in ALL_KEYS if k in values}
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("Настройки сохранены в %s (%d полей)", p, len(out))
