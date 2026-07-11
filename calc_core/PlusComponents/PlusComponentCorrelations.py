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

from calc_core.PlusComponents.CriticalPressure import CriticalPressureCorrelation
from calc_core.PlusComponents.CriticalTemperature import CriticalTemperatureCorrelation
from calc_core.PlusComponents.AcentricFactor import AcentricFactorCorrelation
from calc_core.PlusComponents.CriticalVolume import CriticalVolumeCorrelation
from calc_core.PlusComponents.kWatson import KWatsonCorrelation
from calc_core.PlusComponents.ShiftParameter import ShiftParameterCorrelation
from typing import Dict

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

        self.correlations_config = correlations_config
        self.data = {'molar_mass': M, 'gamma': gamma, 'Tb': Tb}

        if Pc is not None:
            self.data['critical_pressure'] = Pc

        if Tc is not None:
            self.data['critical_temperature'] = Tc

        if Kw is not None:
            self.data['Kw'] = Kw

        if af is not None:
            self.data['acentric_factor'] = af

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
        ValueError
            Если `property_name` неизвестен, метод для него не указан в
            `correlations_config`, метод неизвестен самой корреляции, или
            какой-то требуемый параметр отсутствует в `self.data`.
        """

        if property_name not in property_classes:
            raise ValueError(f"PlusComponentProperties: Unknown property: {property_name}")

        method = self.correlations_config.get(property_name)
        if method is None:
            raise ValueError(f"PlusComponentProperties: No correlation specified for {property_name}")

        calculator = property_classes[property_name]

        try:
            correlation_func = calculator().get_correlation(method.lower())
            required_params = calculator().get_required_params(method.lower())
        except ValueError as e:
            raise ValueError(f"Invalid correlation for {property_name}: {e}")

        # Подготовка параметров
        params = {}
        for param in required_params:
            if param == 'gamma':
                value = self.data.get('gamma')
            elif param == 'Tb_K':
                value = self.data.get('Tb')
            elif param == 'M':
                value = self.data.get('molar_mass')
            elif param == 'Kw':
                value = self.data.get('Kw')
            elif param == 'Tc_K':
                value = self.data.get('critical_temperature')
            elif param == 'Pc_bar':
                value = self.data.get('critical_pressure')
            else:
                value = self.data.get('acentric_factor')

            if value is None:
                raise ValueError(f"Missing required parameter '{param}' for {method}")

            params[param] = value

        return correlation_func(**params)

    def calculate_all(self):
        """
        Считает все 6 свойств псевдокомпонента по очереди (`Kw` →
        `critical_temperature` → `critical_pressure` → `acentric_factor` →
        `critical_volume` → `shift_parameter`) и складывает в `self.data`.
        Порядок важен: некоторые корреляции (например, Kesler-Lee для
        ацентрического фактора) требуют уже посчитанные `Kw`/`Tc`/`Pc`.
        """

        self.data['Kw'] = self.calculate_property('Kw')

        self.data['critical_temperature'] = self.calculate_property('critical_temperature')

        self.data['critical_pressure'] = self.calculate_property('critical_pressure')

        self.data['acentric_factor'] = self.calculate_property('acentric_factor')

        self.data['critical_volume'] = self.calculate_property('critical_volume')

        self.data['shift_parameter'] = self.calculate_property('shift_parameter')
