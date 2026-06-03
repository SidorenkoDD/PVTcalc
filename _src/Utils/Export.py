from abc import ABC, abstractmethod
import sys
import numpy as np
from pathlib import Path
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

root_path = Path(__file__).parent.parent.parent
sys.path.append(str(root_path))

from _src.Utils.Results2 import ResultStore
from _src.CompositionalModel.CompositionalModelV2 import CompositionalModel
# from calculations.Composition.Composition import Composition
# from calculations.CompositionalModel.CompositionalModel import CompositionalModel

# class ExportFacade:
#     def __init__(self):
#         self.E300 = E300()

class Export(ABC):
    @abstractmethod
    def export(self):
        ...

class DB(Export):
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




class E300(Export):
    def __init__(self, model: CompositionalModel):
        self._model = model

    def _write_eos(self, fname):
        fname.write('EOS\n')
        fname.write('--\n')
        fname.write('--Equation of state\n')
        fname.write('--\n')
        if self._model._eos == 'PREOS':
            fname.write('MPR\n')
        elif self._model._eos == 'SRKEOS':
            fname.write('SRK\n')
        fname.write('/\n')


    def _write_n_comps(self, fname):
        fname.write('NCOMPS\n')
        fname.write('--\n')
        fname.write('--Number of components\n')
        fname.write('--\n')
        fname.write(f'{len(self._model._composition._composition.keys())}\n')
        fname.write('/\n')

    def _write_comps(self, fname):
        fname.write('CNAMES\n')
        fname.write('--\n')
        fname.write('--Component Names\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{component}\n')
        fname.write('/\n')

    def _write_zi(self, fname):
        fname.write('ZI\n')
        fname.write('--\n')
        fname.write('--Overall Composition\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{self._model._composition._composition[component]}\n')
        fname.write('/\n')

    def _write_mw(self, fname):
        fname.write('MW\n')
        fname.write('--\n')
        fname.write('--Molecular Weight\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{self._model._composition._composition_data['molar_mass'][component] / 1000}\n')
        fname.write('/\n')

    def _write_tc(self, fname):
        fname.write('TCRIT\n')
        fname.write('--\n')
        fname.write('--Critical Temperature\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{self._model._composition._composition_data['critical_temperature'][component]}\n')
        fname.write('/\n')

    def _write_pc(self, fname):
        fname.write('PCRIT\n')
        fname.write('--\n')
        fname.write('--Critical Pressure\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{self._model._composition._composition_data['critical_pressure'][component] * 10}\n')
        fname.write('/\n')

    def _write_vc(self, fname):
        fname.write('VCRIT\n')
        fname.write('--\n')
        fname.write('--Critical Volume\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{self._model._composition._composition_data['critical_volume'][component]}\n')
        fname.write('/\n')

    def _write_acf(self, fname):
        fname.write('ACF\n')
        fname.write('--\n')
        fname.write('--Acentric factors\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{self._model._composition._composition_data['acentric_factor'][component]}\n')
        fname.write('/\n')

    def _write_bic(self, fname):
        bic_df = self._model._composition.BIPS.iloc[1:]
        mask = np.tril(np.ones_like(bic_df))
        bic_df_diag = bic_df.where(mask.astype(bool))
        fname.write('BIC\n')
        fname.write('--\n')
        fname.write('--Binary interaction coefficients\n')
        fname.write('--\n')
        for row in bic_df_diag.index:
            
            arr = bic_df_diag.loc[row].to_numpy()
            cleared_arr = arr[~np.isnan(arr)]
            result_string = '  '.join(map(str, cleared_arr))

            fname.write(f'{result_string}\n')
        fname.write('/\n')



    def _write_shift(self, fname):
        fname.write('SSHIFT\n')
        fname.write('--\n')
        fname.write('--Eos Volume Shift\n')
        fname.write('--\n')
        for component in list(self._model._composition._composition.keys()):
            fname.write(f'{self._model._composition._composition_data['shift_parameter'][component]}\n')
        fname.write('/\n')

    def _write_info(self, fname):
        fname.write('-- Auto generated file with python\n')


    def export(self, fname:str):
        with open(f'{fname}.inc', 'w') as result_file:
            result_file.write('--E300.inc\n')
            self._write_info(result_file)
            self._write_eos(result_file)
            self._write_n_comps(result_file)
            self._write_comps(result_file)
            self._write_mw(result_file)
            self._write_tc(result_file)
            self._write_pc(result_file)
            self._write_vc(result_file)
            self._write_acf(result_file)
            self._write_shift(result_file)
            self._write_zi(result_file)
            self._write_bic(result_file)

