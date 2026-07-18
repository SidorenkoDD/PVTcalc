"""
История расчётов в памяти (`ResultStore`) — не путать с `Utils/Export.py`
(персистентность моделей в JSON). Используется `CompositionalModel`, куда
`ResultStore` попадает как `result_store_object` и логирует каждый вызов
`.flash(P, T)`.
"""

from dataclasses import dataclass, field
from typing import Any, Dict
from datetime import datetime
from typing import List, Optional
from pathlib import Path
import uuid
import pickle





@dataclass
class CalcResult:
    """Один залогированный расчёт: модуль-источник, входные параметры, результат, метка времени, уникальный id."""
    module : str
    params : Dict[str, Any]
    data : Any
    timestamp : datetime = field(default_factory=datetime.now)
    id : str = field(default_factory= lambda: str(uuid.uuid4()))

class ResultStore:
    """История расчётов в оперативной памяти (список `CalcResult`), с поиском по id/модулю."""

    def __init__(self):
        self._results: List[CalcResult] = []

    def add(self, module: str, params: Dict, data: Any) -> str:
        """
        Добавляет запись в историю.

        Parameters
        ----------
        module : str
            Имя модуля-источника (например, `'Flash'`).
        params : Dict
            Входные параметры расчёта (например, `{'P': P, 'T': T}`).
        data : Any
            Результат расчёта — как есть, без сериализации.

        Returns
        -------
        str
            Сгенерированный `id` новой записи (uuid4).
        """
        res = CalcResult(module=module,
                         params=params,
                         data=data)
        self._results.append(res)
        return res.id
    
    def get(self, result_id: str) -> Optional[CalcResult]:
        """Находит запись по `id`, либо `None`, если такой нет."""
        return next((r for r in self._results if r.id == result_id), None)

    def get_by_module( self, module: str) -> List[CalcResult]:
        """Все записи для заданного модуля-источника (в порядке добавления)."""
        return [r for r in self._results if r.module == module]

    def save(self, fpath:str):
        """
        Сохраняет всю историю в файл через `pickle`.

        **Известный баг** (CLAUDE.md Known Issues): файл открывается в
        режиме `'rb'` (чтение) вместо `'wb'` (запись) — вызов гарантированно
        падает с `io.UnsupportedOperation` при первом же обращении. Не
        исправлено намеренно, не чинить попутно без согласования.

        Parameters
        ----------
        fpath : str
        """
        Path(fpath).parent.mkdir(parents=True, exist_ok=True)
        with open(fpath, 'wb') as f:
            pickle.dump(self._results, f)
