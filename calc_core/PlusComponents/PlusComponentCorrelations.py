"""
Агрегатор корреляций свойств псевдокомпонентов C7+.

Диспетчеризует расчёт по строковым ключам метода (`correlations_config`,
например `{'critical_temperature': 'pedersen', ...}`) к соответствующему
классу-корреляции в этой же папке (`CriticalTemperature.py`,
`CriticalPressure.py`, `AcentricFactor.py`, `CriticalVolume.py`,
`kWatson.py`, `ShiftParameter.py` — по одному файлу на свойство). Каждый
такой класс сам знает, какие параметры (`gamma`, `M`, `Tb_K`, `Kw`, `Tc_K`,
`Pc_bar`, `af`) нужны конкретному методу (`get_required_params`) — эта
подготовка параметров сделана в `calculate_property`. Используется из
`Composition.evaluate_composition_data` для компонент с `c7_plus_flag=True`.
"""

import math
from typing import Dict

from calc_core.PlusComponents.AcentricFactor import AcentricFactorCorrelation
from calc_core.PlusComponents.CriticalPressure import CriticalPressureCorrelation
from calc_core.PlusComponents.CriticalTemperature import CriticalTemperatureCorrelation
from calc_core.PlusComponents.CriticalVolume import CriticalVolumeCorrelation
from calc_core.PlusComponents.kWatson import KWatsonCorrelation
from calc_core.PlusComponents.ShiftParameter import ShiftParameterCorrelation
from calc_core.Utils.Errors import InputValidationError

property_classes = {
    'critical_temperature': CriticalTemperatureCorrelation,
    'critical_pressure': CriticalPressureCorrelation,
    'acentric_factor': AcentricFactorCorrelation,
    'critical_volume': CriticalVolumeCorrelation,
    'Kw': KWatsonCorrelation,
    'shift_parameter': ShiftParameterCorrelation
}

class PlusComponentProperties:
    """
    Main class to calculate C7+ components properties.
    """

    _POSITIVE_PROPERTIES = {
        'critical_temperature',
        'critical_pressure',
        'critical_volume',
        'Kw',
    }

    @staticmethod
    def _finite_number(value, name: str, *, positive: bool = False) -> float:
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise InputValidationError(
                f'{name} должно быть числом. Передано: {value!r}'
            ) from exc
        if isinstance(value, bool) or not math.isfinite(number):
            raise InputValidationError(
                f'{name} должно быть конечным числом. Передано: {value!r}'
            )
        if positive and number <= 0.0:
            raise InputValidationError(
                f'{name} должно быть больше 0. Передано: {value!r}'
            )
        return number

    def __init__(self, M: float, gamma: float, Tb: float, correlations_config: Dict[str, str], *, Tc: float = None,
                 Pc: float = None, Kw: float = None, af: float = None):
        """
        Parameters
        ----------
        M : float
            Молярная масса псевдокомпонента, г/моль.
        gamma : float
            Относительная плотность.
        Tb : float
            Температура кипения, K.
        correlations_config : Dict[str, str]
            `{имя_свойства: имя_метода}`, например
            `{'critical_temperature': 'pedersen', 'critical_pressure': 'rizari_daubert',
            'acentric_factor': 'riazi al sahhaf', 'critical_volume': 'hall yarborough',
            'Kw': 'k watson', 'shift_parameter': 'jhaveri youngren'}`.
            Ключи — из `property_classes`, значения — методы, известные
            соответствующему классу-корреляции (см. `get_correlation`
            каждого из них за списком допустимых значений).
        Tc, Pc, Kw, af : float, optional
            Заранее известные значения (если какое-то свойство уже
            посчитано извне и пересчитывать через корреляцию не нужно) —
            подставляются в `self.data` напрямую, минуя корреляции.
        """

        if not isinstance(correlations_config, dict):
            raise InputValidationError('correlations_config должен быть словарём.')

        self.correlations_config = correlations_config
        self.data = {
            'molar_mass': self._finite_number(M, 'M', positive=True),
            'gamma': self._finite_number(gamma, 'gamma', positive=True),
            'Tb': self._finite_number(Tb, 'Tb', positive=True),
        }

        if Pc is not None:
            self.data['critical_pressure'] = self._finite_number(Pc, 'Pc', positive=True)

        if Tc is not None:
            self.data['critical_temperature'] = self._finite_number(Tc, 'Tc', positive=True)

        if Kw is not None:
            self.data['Kw'] = self._finite_number(Kw, 'Kw', positive=True)

        if af is not None:
            self.data['acentric_factor'] = self._finite_number(af, 'af')

    def calculate_property(self, property_name: str) -> float:
        """
        Универсальный метод расчета любого свойства.

        Смотрит выбранный метод в `self.correlations_config[property_name]`,
        находит нужную функцию-корреляцию (`property_classes[property_name].get_correlation(method)`),
        собирает требуемые ей параметры из `self.data` и вызывает.

        Parameters
        ----------
        property_name : str
            Один из ключей `property_classes` (`'critical_temperature'`,
            `'critical_pressure'`, `'acentric_factor'`, `'critical_volume'`,
            `'Kw'`, `'shift_parameter'`).

        Returns
        -------
        float

        Raises
        ------
        ValueError or InputValidationError
            Если `property_name` неизвестен, метод для него не указан в
            `correlations_config`, метод неизвестен самой корреляции, или
            какой-то требуемый параметр отсутствует/нечисловой в `self.data`.
        """

        if property_name not in property_classes:
            raise ValueError(f"PlusComponentProperties: Unknown property: {property_name}")

        method = self.correlations_config.get(property_name)
        if method is None:
            raise ValueError(f"PlusComponentProperties: No correlation specified for {property_name}")
        if not isinstance(method, str):
            raise InputValidationError(
                f"Correlation method for {property_name} must be a string. Got: {method!r}"
            )

        calculator = property_classes[property_name]

        try:
            correlation_func = calculator().get_correlation(method.lower())
            required_params = calculator().get_required_params(method.lower())
        except ValueError as e:
            raise ValueError(f"Invalid correlation for {property_name}: {e}")

        # Подготовка параметров
        param_to_data_key = {
            'gamma': 'gamma',
            'Tb_K': 'Tb',
            'M': 'molar_mass',
            'Kw': 'Kw',
            'Tc_K': 'critical_temperature',
            'Pc_bar': 'critical_pressure',
            'af': 'acentric_factor',
        }
        params = {}
        for param in required_params:
            data_key = param_to_data_key.get(param)
            if data_key is None:
                raise ValueError(
                    f"Correlation '{method}' declares unsupported parameter '{param}'"
                )
            value = self.data.get(data_key)

            if value is None:
                raise ValueError(f"Missing required parameter '{param}' for {method}")

            params[param] = value

        result = self._finite_number(
            correlation_func(**params),
            f'{property_name} ({method})',
            positive=property_name in self._POSITIVE_PROPERTIES,
        )
        return result

    def calculate_all(self):
        """
        Заполняет до 6 свойств псевдокомпонента по очереди (`Kw` →
        `critical_temperature` → `critical_pressure` → `acentric_factor` →
        `critical_volume` → `shift_parameter`) и складывает в `self.data`.
        Порядок важен: некоторые корреляции (например, Kesler-Lee для
        ацентрического фактора) требуют уже посчитанные `Kw`/`Tc`/`Pc`.
        Значения `Tc`/`Pc`/`Kw`/`af`, переданные явно в `__init__`, считаются
        исходными данными и пропускаются без пересчёта.
        """

        calculation_order = (
            'Kw',
            'critical_temperature',
            'critical_pressure',
            'acentric_factor',
            'critical_volume',
            'shift_parameter',
        )
        for property_name in calculation_order:
            # Явно переданные Tc/Pc/Kw/af являются входными данными и не
            # должны перезаписываться корреляцией (см. контракт __init__).
            if property_name not in self.data:
                self.data[property_name] = self.calculate_property(property_name)
