"""
Персистентный слой: сохранение/загрузка "снэпшотов" моделей флюидов.

Пишет в JSON-файл (обычно `models.json` в корне репозитория) — состав,
`composition_data`, тип EOS и историю результатов расчётов для именованных
моделей. Читается обратно через `Composition.from_db(path)`. **Важно**: это
не независимый эталонный датасет (например, PVTSim) — сюда пишутся
результаты расчётов самого движка, см. CLAUDE.md.
"""

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional
# from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
# from calc_core.Utils.Results import ResultStore

logger = logging.getLogger(__name__)


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
        
        # ВАЖНО: Сразу загружаем существующие данные при инициализации,
        # чтобы новые записи добавлялись к старым, а не заменяли их.
        self.load()

    def export(self, model_id: str,
               name: str,
               composition: Dict[str, Any],
               composition_data: Dict[str, Any],
               eos: Any,
               results: Optional[Any] = None, # Сделали Optional, так как может быть None
               field: str = None):
        """Кладёт модель в оперативную память 'базы'."""
        
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
        self._db[model_id] = {
            "Model_name": name,
            "Field": field,
            "composition": composition,
            "composition_data": composition_data,
            "eos": eos, 
            "created_at": datetime.now().isoformat(),
            "results": res_serialized
        }

    def save(self):
        """Записывает все накопленные модели в файл."""
        self.filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(self.filepath, "w", encoding="utf-8") as f:
            json.dump(self._db, f, indent=2, default=str, ensure_ascii=False)

    def load(self):
        """Загружает файл в текущее содержимое. Если файла нет, оставляет self._db пустым."""
        if self.filepath.exists():
            try:
                with open(self.filepath, "r", encoding="utf-8") as f:
                    self._db = json.load(f)
            except json.JSONDecodeError:
                # На случай, если файл поврежден, начинаем с чистого листа
                logger.warning("%s поврежден. Начинаем с пустой базы.", self.filepath)
                self._db = {}

    @staticmethod
    def _safe_json(data: Any) -> Any:
        """Превращает numpy/pandas в list/dict, чтобы JSON не падал."""
        if hasattr(data, "tolist"): return data.tolist()
        if hasattr(data, "to_dict"): return data.to_dict(orient="list")
        return data

    # --- ДОПОЛНИТЕЛЬНО: Удобный метод "всё в одном" ---
    def export_and_save(self, model_id: str, name: str, composition: Dict[str, Any],
                        composition_data: Dict[str, Any], eos: Any, 
                        results: Optional[Any] = None, field: str = None):
        """Добавляет модель в память и сразу сохраняет на диск."""
        self.export(model_id, name, composition, composition_data, eos, results, field)
        self.save()