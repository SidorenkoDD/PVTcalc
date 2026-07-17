"""
Экспорт активной модели в внешние форматы (для GUI).

Пока поддержан один формат — Eclipse 300 (`.inc`) через
`calc_core.Utils.E300Export.E300Exporter`. Список форматов (`FORMATS`)
описывает подпись, расширение и опции формата (для окна экспорта). Не
импортирует DearPyGui.
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.Utils.E300Export import E300Exporter

logger = logging.getLogger(__name__)

# Описание форматов экспорта для окна: подпись, расширение, опции (EOS и т.п.).
FORMATS: dict[str, dict] = {
    "e300": {
        "label": "Eclipse 300 (.inc)",
        "ext": ".inc",
        "eos_choices": ["MPR", "SRK"],  # Брусиловский -> выбор EOS для дека
    },
}


def format_labels() -> list[str]:
    """Подписи форматов для выпадающего списка (в порядке словаря)."""
    return [meta["label"] for meta in FORMATS.values()]


def format_key_by_label(label: str) -> str:
    """Ключ формата по его подписи (обратное к `format_labels`)."""
    for key, meta in FORMATS.items():
        if meta["label"] == label:
            return key
    return next(iter(FORMATS))


def export_model(composition: Composition, fmt: str, path: str,
                 eos_keyword: str = "MPR") -> str:
    """
    Экспортирует состав в формат `fmt` в файл `path`. Возвращает путь.

    Parameters
    ----------
    fmt : str
        Ключ формата из `FORMATS` (пока только `"e300"`).
    eos_keyword : str
        Для E300 — уравнение состояния дека (`MPR`/`SRK`).
    """
    if fmt == "e300":
        return E300Exporter(composition, eos_keyword=eos_keyword).export(path)
    raise ValueError(f"Неизвестный формат экспорта: {fmt}")
