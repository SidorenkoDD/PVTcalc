"""
Заготовка под чтение сохранённых моделей "в обратную сторону" от `Export.py`.

`Export2.py::ModelJSONDB` пишет `Composition` → `models.json`; по замыслу
`DBModelImport` должен делать обратное — читать `models.json` и строить
`Composition`/список моделей программно (сейчас эта функциональность
фактически покрыта `Composition.from_db(...)`, который встроен прямо в
Composition.py). Не реализовано — класс с пустым `__init__`, нигде не
вызывается. Автор подтвердил, что модуль понадобится в будущем.
"""

from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any
from calc_core.Utils.Results import ResultStore


class DBModelImport:
    """Не реализовано (заготовка) — см. module docstring."""

    def __init__(self):
        pass