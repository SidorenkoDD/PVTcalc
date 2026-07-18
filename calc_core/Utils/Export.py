"""
Персистентный слой: сохранение/загрузка "снэпшотов" моделей флюидов.

Пишет в JSON-файл (обычно `models.json` в корне репозитория) — состав,
`composition_data`, тип EOS и историю результатов расчётов для именованных
моделей. Читается обратно через `Composition.from_db(path)`. **Важно**: это
не независимый эталонный датасет (например, PVTSim) — сюда пишутся
результаты расчётов самого движка, см. CLAUDE.md.
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from calc_core.Utils.ModelStore import (
    StoreRevision,
    read_model_store,
    store_revision,
    write_model_store,
)

# from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
# from calc_core.Utils.Results import ResultStore

class ModelJSONDB:
    """Простая JSON-"база" в памяти (`self._db`), синхронизируемая с файлом на диске через `load()`/`save()`."""

    def __init__(self, filepath: str = "models.json"):
        """
        Parameters
        ----------
        filepath : str, optional
            Путь к JSON-файлу. По умолчанию `"models.json"`. Существующий
            файл сразу загружается в `self._db` (см. `load()`), чтобы
            `export()` добавлял модели к уже сохранённым, а не затирал их.
        """
        self.filepath = Path(filepath)
        self._db: Dict[str, Dict[str, Any]] = {}
        self._source_revision: StoreRevision = None
        
        # ВАЖНО: Сразу загружаем существующие данные при инициализации,
        # чтобы новые записи добавлялись к старым, а не заменяли их.
        self.load()

    def export(self, model_id: str,
               name: str,
               composition: Dict[str, Any],
               composition_data: Dict[str, Any],
               eos: Any,
               results: Optional[Any] = None, # Сделали Optional, так как может быть None
               field: Optional[str] = None,
               t_res: Optional[float] = None,
               correlations: Optional[Dict[str, Any]] = None,
               project_id: Optional[str] = None,
               project_name: Optional[str] = None):
        """Кладёт модель в оперативную память 'базы'.

        `t_res` — пластовая температура, K (опционально). Если задана —
        пишется в запись ключом `"T_res"`; старые записи без ключа остаются
        валидными (`Composition.from_db` для них использует default_T).
        """

        # ЗАЩИТА: Если results == None, используем пустой список, чтобы не было ошибки в цикле
        safe_results = results if results is not None else []

        res_serialized = [
            {
                "id": r.id,
                "module": r.module,
                "params": r.params,
                "data": self._safe_json(r.data),
                "timestamp": r.timestamp.isoformat(),
            }
            for r in safe_results
        ]

        # Если model_id уже существует, он будет обновлен (перезаписан).
        # Если нет — будет создан новый.
        previous = self._db.get(model_id, {})
        now = datetime.now().isoformat()
        record = {
            "Model_name": name,
            "Field": field,
            "composition": composition,
            "composition_data": composition_data,
            "eos": eos,
            "created_at": previous.get("created_at", now),
            "updated_at": now,
            "results": res_serialized,
        }
        if t_res is not None:
            record["T_res"] = float(t_res)
        if correlations is not None:
            record["correlations"] = correlations
        if project_id is not None:
            record["project_id"] = project_id
        if project_name is not None:
            record["project_name"] = project_name
        self._db[model_id] = record

    def save(self):
        """Валидирует и атомарно пишет store с lock/backup/concurrency check."""
        write_model_store(
            self.filepath,
            self._db,
            expected_revision=self._source_revision,
        )
        self._source_revision = store_revision(self.filepath)

    def load(self):
        """Загружает файл в текущее содержимое. Если файла нет, оставляет self._db пустым."""
        self._db = read_model_store(self.filepath, missing_ok=True)
        self._source_revision = store_revision(self.filepath)

    @staticmethod
    def _safe_json(data: Any) -> Any:
        """Превращает numpy/pandas в list/dict, чтобы JSON не падал."""
        if hasattr(data, "tolist"):
            return data.tolist()
        if hasattr(data, "to_dict"):
            return data.to_dict(orient="list")
        return data

    # --- ДОПОЛНИТЕЛЬНО: Удобный метод "всё в одном" ---
    def export_and_save(self, model_id: str, name: str, composition: Dict[str, Any],
                        composition_data: Dict[str, Any], eos: Any,
                        results: Optional[Any] = None, field: Optional[str] = None,
                        t_res: Optional[float] = None,
                        correlations: Optional[Dict[str, Any]] = None,
                        project_id: Optional[str] = None,
                        project_name: Optional[str] = None):
        """Добавляет модель в память и сразу сохраняет на диск."""
        self.export(model_id, name, composition, composition_data, eos, results,
                    field, t_res, correlations, project_id, project_name)
        self.save()
