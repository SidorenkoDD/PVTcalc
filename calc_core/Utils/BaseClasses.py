"""
Набор абстрактных заготовок, из которых в живом коде используется только `Reader`.

`Reader` — родитель `Utils/JsonDBReader.py::JsonDBReader` (реально работает).
`ModuleProperty`, `PhaseStabilityTest`, `Calculator`, `CalculationModule` —
нигде не наследуются и не инстанцируются в текущей кодовой базе; похоже на
незавершённый набросок общей архитектуры "модуль с зарегистрированными
калькуляторами", от которого в итоге остался только `Reader`.
"""

from abc import abstractmethod, ABC
#from calculations.Composition.Composition import Composition
from typing import Dict, TypeVar, Generic
from calc_core.EOS.BaseEOS import EOS


T = TypeVar('T')

class ModuleProperty(Generic[T]):
    """Дескриптор для модулей с type hints. Нигде не используется (см. module docstring)."""
    def __get__(self, obj, owner) -> T:
        if obj is None:
            return self
        return obj._modules[self._name]

    def __set_name__(self, owner, name):
        self._name = name




class PhaseStabilityTest:
    '''
    Абстрактный класс для stability test. Нигде не наследуется — живой
    `PhaseStability.TwoPhaseStabilityTest` от него не зависит.
    '''
    def __init__(self, composition, p, t, eos: EOS | str):
        self.composition = composition
        self.p = p
        self.t = t
        self.eos = eos


    @abstractmethod
    def calculate_phase_stability(self) -> dict:
        pass



class Calculator:
    '''
    Абстрактный класс для расчетного подмодуля (например, TwoPhaseFlash).
    Нигде не наследуется.
    '''

    def __init__(self, composition, eos):
        self.composition = composition
        self.eos = eos
        self.last_result = None

    def calculate(self, conditions):
        raise NotImplementedError


class CalculationModule:
    '''
    Абстрактный класс для расчетного модуля (например, Flash) с реестром
    именованных калькуляторов, доступных через `__getattr__`. Нигде не наследуется.
    '''

    def __init__(self):
        self._calculators = {}

    def register_calculator(self, name, calculator):
        """Регистрирует калькулятор под именем `name`, делая его доступным как атрибут (см. `__getattr__`)."""
        self._calculators[name] = calculator

    def __getattr__(self, name):
        """Возвращает зарегистрированный калькулятор по имени; иначе `AttributeError`."""
        if name in self._calculators:
            return self._calculators[name]
        raise AttributeError(f"No calculator '{name}' in module")

class Reader(ABC):
    '''
    Abstract class for data reader — единственный живой класс в этом файле.
    Наследуется `Utils/JsonDBReader.py::JsonDBReader`.
    '''
    @abstractmethod
    def load_database(self):
        pass