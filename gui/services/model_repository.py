"""
Репозиторий сохранённых моделей флюидов.

Единственная точка контакта GUI с персистентным стором `models.json`.
Список моделей и загрузка составов идут через уже существующий
`Composition.from_db(...)` (см. calc_core/Composition/Composition.py) — он
отдаёт прокси с `list_models()` и построением `Composition` по ключу без
пересчёта корреляций (используются сохранённые числа). Метаданные для
списка (месторождение, EOS, число компонентов) читаются напрямую из JSON,
т.к. прокси их не выставляет.

Заготовка `calc_core/Utils/Import.py::DBModelImport` для чтения намеренно
не используется — она пустая, а `from_db` уже покрывает нужное (см.
docs/GUI.md, заметка 2026-07-16).
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from calc_core.Composition.Composition import Composition

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ModelSummary:
    """Лёгкая сводка модели для дерева/списка — без загрузки состава."""

    model_id: str
    title: str
    field_name: Optional[str]
    eos: Optional[str]
    n_components: int


class ModelRepository:
    """Чтение списка и составов моделей из `models.json`."""

    def __init__(self, db_path: str = "models.json", t_res: float = 373.15):
        """
        Parameters
        ----------
        db_path : str
            Путь к файлу снэпшотов моделей. По умолчанию `models.json`.
        t_res : float
            Пластовая температура, проставляемая загруженным составам
            (в файле не хранится — см. `Composition.from_db`).
        """
        self._db_path = Path(db_path)
        self._t_res = t_res

    @property
    def db_path(self) -> Path:
        return self._db_path

    def list_models(self) -> list[ModelSummary]:
        """
        Возвращает сводки всех моделей в файле. Пустой список, если файла
        нет (в отличие от `from_db`, который бросает FileNotFoundError) —
        удобнее для UI: показать пустое дерево, а не падать.
        """
        if not self._db_path.exists():
            logger.warning("Файл моделей не найден: %s", self._db_path)
            return []

        with open(self._db_path, "r", encoding="utf-8") as f:
            raw = json.load(f)

        summaries: list[ModelSummary] = []
        for model_id, rec in raw.items():
            summaries.append(
                ModelSummary(
                    model_id=model_id,
                    title=rec.get("Model_name") or model_id,
                    field_name=rec.get("Field"),
                    eos=str(rec.get("eos")) if rec.get("eos") is not None else None,
                    n_components=len(rec.get("composition", {})),
                )
            )
        return summaries

    def load_composition(self, model_id: str) -> Composition:
        """Строит `Composition` для модели через `Composition.from_db`."""
        db = Composition.from_db(str(self._db_path), T_res=self._t_res)
        if model_id not in db.list_models():
            raise KeyError(
                f"Модель '{model_id}' не найдена в {self._db_path}. "
                f"Доступные: {db.list_models()}"
            )
        composition = getattr(db, model_id)
        logger.info("Состав модели '%s' загружен", model_id)
        return composition
