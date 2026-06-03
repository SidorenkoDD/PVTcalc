from _src.CompositionalModel.CompositionalModelV2 import CompositionalModel
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any
from _src.Utils.Results2 import ResultStore  # твой класс

class ModelJSONDB:
    def __init__(self, filepath: str = "models.json"):
        self.filepath = Path(filepath)
        self._db: Dict[str, Dict[str, Any]] = {}  # id_модели → её данные

    def export(self, model_id: str,
               name: str,
               composition: Dict[str, Any],
               composition_data : Dict[str, Any],
               eos,
               results: ResultStore = None,
               field: str = None):
        

        """Кладёт модель в оперативную память 'базы'."""
        # Превращаем результаты в обычные словари/списки
        res_serialized = [
            {
                "id": r.id,
                "module": r.module,
                "params": r.params,
                "data": self._safe_json(r.data),
                "timestamp": r.timestamp.isoformat(),
            }

            for r in results
        ]

        self._db[model_id] = {
            "Model_name": name,
            "Field" : field,
            "composition": composition,
            "composition_data" : composition_data,
            "eos" : eos, 
            "created_at": datetime.now().isoformat(),
            "results": res_serialized
        }
    
    
    def save(self):
        """Записывает все накопленные модели в один файл."""
        self.filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(self.filepath, "w", encoding="utf-8") as f:
            # indent=2 для читаемости, default=str страховка от редких не-сериализуемых объектов
            json.dump(self._db, f, indent=2, default=str, ensure_ascii=False)

    def load(self):
        """Загружает файл вместо текущего содержимого."""
        if self.filepath.exists():
            with open(self.filepath, "r", encoding="utf-8") as f:
                self._db = json.load(f)

    @staticmethod
    def _safe_json(data: Any) -> Any:
        """Превращает numpy/pandas в list/dict, чтобы JSON не падал."""
        if hasattr(data, "tolist"): return data.tolist()
        if hasattr(data, "to_dict"): return data.to_dict(orient="list")
        return data