'''
module to import composition and create dict

Альтернатива `Composition.from_db(...)` — загрузка состава из внешнего
Excel-файла (двухколоночный лист `key`/`value`), а не из `models.json`.
Используется в `test_notebook.ipynb`, раздел "Старая загрузка составов",
с жёстко прописанным путём к файлу на конкретной машине автора — не
портируемо между машинами как есть.
'''
from abc import ABC, abstractmethod
from calc_core.Utils.Errors import InvalidExcelComponentType, InvalidExcelValueType
from calc_core.Composition.Composition import Composition
from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
import pandas as pd

class CompositionLoader(ABC):
    """Абстрактный загрузчик состава из внешнего файла."""

    def __init__(self, filepath):
        """
        Parameters
        ----------
        filepath : str
            Путь к файлу-источнику.
        """
        self._filepath = filepath
    @abstractmethod
    def load(self):
        pass



class CompositionExcelLoader(CompositionLoader):
    """Загрузка состава из Excel-файла с двухколоночным листом `key`/`value`."""

    def __init__(self, filepath):
        super().__init__(filepath)

    def load(self, header: bool, sheet : str):
        """
        Читает лист Excel и превращает его в `{имя_компонента: мольная_доля}`.

        Parameters
        ----------
        header : bool
            `True` — лист читается со строкой заголовков
            (`pd.read_excel(..., sheet_name=sheet)`); `False` — без
            заголовков (`header=None`), но при этом параметр `sheet`
            игнорируется (не передаётся в вызов `read_excel`) — читается
            активный/первый лист файла независимо от значения `sheet`.
        sheet : str
            Имя листа — используется только при `header=True`.

        Returns
        -------
        dict
            `{имя_компонента: мольная_доля}`.

        Raises
        ------
        InvalidExcelValueType
            Если колонка значений не числовая.
        InvalidExcelComponentType
            Если колонка ключей не строковая.
        """
        if header:
            df = pd.read_excel(io = self._filepath, sheet_name= sheet)
            df.columns = ["key", "value"]
            if not pd.api.types.is_numeric_dtype(df['value']):
                raise InvalidExcelValueType('Value not numeric in Excel!')
            if not pd.api.types.is_string_dtype(df['key']):
                raise InvalidExcelComponentType('Value not string in Excel!')
            dict_from_df = df.set_index('key')['value'].to_dict()
        else:
            df = pd.read_excel(io = self._filepath,
                            header = None)
            df.columns = ["key", "value"]
            if not pd.api.types.is_numeric_dtype(df['value']):
                raise InvalidExcelValueType('Value not numeric in Excel!')
            if not pd.api.types.is_string_dtype(df['key']):
                raise InvalidExcelComponentType('Value not string in Excel!')
            dict_from_df = df.set_index('key')['value'].to_dict()

        return dict_from_df