from dataclasses import dataclass, field
from typing import Any, Dict
from datetime import datetime
from typing import List, Optional
from pathlib import Path
import uuid
import pickle



@dataclass
class CalcResult:
    module : str
    params : Dict[str, Any]
    data : Any
    timestamp : datetime = field(default_factory=datetime.now)
    id : str = field(default_factory= lambda: str(uuid.uuid4()))

class ResultStore:
    def __init__(self):
        self._results: List[CalcResult] = []

    def add(self, module: str, params: Dict, data: Any) -> str:
        res = CalcResult(module=module,
                         params=params,
                         data=data)
        self._results.append(res)
        return res.id
    
    def get(self, result_id: str) -> Optional[CalcResult]:
        return next((r for r in self._results if r.id == result_id), None)
    
    def get_by_module( self, module: str) -> List[CalcResult]:
        return [r for r in self._results if r.module == module]
    
    def save(self, fpath:str):
        Path(fpath).parent.mkdir(parents=True, exist_ok=True)
        with open(fpath, 'rb') as f:
            pickle.dump(self._results, f)
