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

import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from calc_core.Composition.Composition import Composition
from calc_core.Utils.ModelStore import read_model_store, update_model_store

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ModelSummary:
    """Лёгкая сводка модели для дерева/списка — без загрузки состава."""

    model_id: str
    title: str
    field_name: Optional[str]
    eos: Optional[str]
    n_components: int
    created_at: Optional[str] = None
    t_res: Optional[float] = None
    project_id: Optional[str] = None
    project_name: Optional[str] = None
    # лёгкий срез results из models.json (без data): [{"module", "timestamp"}]
    results_brief: tuple = ()


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
            без собственного ``T_res`` (для legacy-записей).
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

        raw = read_model_store(self._db_path)

        summaries: list[ModelSummary] = []
        for model_id, rec in raw.items():
            results_brief = tuple(
                {"module": r.get("module"), "timestamp": r.get("timestamp")}
                for r in rec.get("results", []) if isinstance(r, dict)
            )
            t_res = rec.get("T_res")
            project_id = rec.get("project_id")
            if not isinstance(project_id, str) or not project_id.strip():
                project_id = model_id
            else:
                project_id = project_id.strip()
            project_name = rec.get("project_name")
            if not isinstance(project_name, str) or not project_name.strip():
                project_name = rec.get("Model_name") or model_id
            summaries.append(
                ModelSummary(
                    model_id=model_id,
                    title=rec.get("Model_name") or model_id,
                    field_name=rec.get("Field"),
                    eos=str(rec.get("eos")) if rec.get("eos") is not None else None,
                    n_components=len(rec.get("composition", {})),
                    created_at=rec.get("created_at"),
                    t_res=float(t_res) if t_res is not None else None,
                    project_id=project_id,
                    project_name=project_name,
                    results_brief=results_brief,
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

    def load_correlations(self, model_id: str) -> Optional[dict]:
        """Возвращает сохранённый выбор корреляций C7+ модели, если он есть."""
        if not self._db_path.exists():
            return None
        raw = read_model_store(self._db_path)
        rec = raw.get(model_id)
        if not isinstance(rec, dict):
            return None
        value = rec.get("correlations")
        return dict(value) if isinstance(value, dict) else None

    def save_composition(self, model_id: str, composition: Composition,
                         correlations: Optional[dict] = None) -> None:
        """Атомарно сохраняет отредактированный состав, сохраняя метаданные.

        Перед записью ротируются ``.bak``, ``.bak.1`` и ``.bak.2``.
        ``created_at`` и история результатов не сбрасываются; меняется только
        редактируемая часть модели и ``updated_at``.
        """
        def update(raw) -> None:
            if model_id not in raw:
                raise KeyError(f"Модель '{model_id}' не найдена в {self._db_path}")
            rec = raw[model_id]
            rec["composition"] = dict(composition.composition)
            rec["composition_data"] = composition.composition_data
            rec["eos"] = composition.eos_name.value
            rec["T_res"] = float(composition.T)
            rec["correlations"] = dict(correlations or {})
            rec["updated_at"] = datetime.now().isoformat(timespec="seconds")

        update_model_store(self._db_path, update)
        logger.info("Модель '%s' сохранена в %s", model_id, self._db_path)
