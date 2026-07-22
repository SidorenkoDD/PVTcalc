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

import hashlib
import json
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

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
    # Время последнего полного сохранения графа/результатов в records["workspace"].
    workspace_saved_at: Optional[str] = None


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
                dict(r)
                for r in rec.get("results", []) if isinstance(r, dict)
            )
            workspace = rec.get("workspace")
            workspace_saved_at = (workspace.get("saved_at")
                                  if isinstance(workspace, dict) else None)
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
                    workspace_saved_at=(workspace_saved_at
                                        if isinstance(workspace_saved_at, str)
                                        else None),
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

    def load_workspace(self, model_id: str) -> Optional[dict]:
        """Возвращает сохранённый GUI-workspace модели, если он есть.

        Результаты являются частью записи модели, а не `gui_session.json`.
        Обёртка `workspace` оставляет место для версии, provenance и будущих
        миграций, тогда как наружу отдаётся только снимок codec-а.
        """
        if not self._db_path.exists():
            return None
        rec = read_model_store(self._db_path).get(model_id)
        if not isinstance(rec, dict):
            return None
        stored = rec.get("workspace")
        if not isinstance(stored, dict):
            return None
        snapshot = stored.get("snapshot")
        return dict(snapshot) if isinstance(snapshot, dict) else None

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

    def save_model_snapshot(self, model_id: str, composition: Composition,
                            correlations: Optional[dict], workspace: dict) -> None:
        """Атомарно сохраняет состав, граф workspace и индекс его результатов.

        `workspace` уже сериализован `workspace_codec`; исходные таблицы и
        результаты остаются в нём целиком, а `results` служит компактной
        сводкой для Projects без загрузки тяжёлых данных.
        """
        if not isinstance(workspace, dict):
            raise TypeError("workspace must be a JSON object")
        now = datetime.now().isoformat(timespec="seconds")
        composition_data = composition.composition_data
        provenance_payload: dict[str, Any] = {
            "composition": dict(composition.composition),
            "composition_data": composition_data,
            "eos": composition.eos_name.value,
            "T_res": float(composition.T),
            "correlations": dict(correlations or {}),
        }
        provenance = hashlib.sha256(
            json.dumps(provenance_payload, sort_keys=True,
                       ensure_ascii=False, allow_nan=False).encode("utf-8")
        ).hexdigest()

        results: list[dict[str, Any]] = []
        raw_nodes = workspace.get("nodes")
        if isinstance(raw_nodes, list):
            for raw_node in raw_nodes:
                if not isinstance(raw_node, dict):
                    continue
                kind = raw_node.get("kind")
                if kind == "composition" or not isinstance(kind, str):
                    continue
                params = raw_node.get("params")
                results.append({
                    "node_id": raw_node.get("node_id"),
                    "module": raw_node.get("title") or kind,
                    "kind": kind,
                    "status": raw_node.get("status"),
                    "has_result": raw_node.get("result") is not None,
                    "experiment_kind": (params.get("kind") if isinstance(params, dict)
                                        and isinstance(params.get("kind"), str) else None),
                    "timestamp": now,
                    "provenance": provenance,
                })

        def update(raw) -> None:
            if model_id not in raw:
                raise KeyError(f"Модель '{model_id}' не найдена в {self._db_path}")
            rec = raw[model_id]
            rec["composition"] = dict(composition.composition)
            rec["composition_data"] = composition_data
            rec["eos"] = composition.eos_name.value
            rec["T_res"] = float(composition.T)
            rec["correlations"] = dict(correlations or {})
            rec["workspace"] = {
                "schema_version": 1,
                "codec_version": workspace.get("schema_version"),
                "saved_at": now,
                "algorithm": "PVTcalc GUI workspace",
                "provenance": {"composition_sha256": provenance},
                "snapshot": workspace,
            }
            rec["results"] = results
            rec["updated_at"] = now

        update_model_store(self._db_path, update)
        logger.info("Модель '%s' и её workspace сохранены в %s", model_id,
                    self._db_path)
