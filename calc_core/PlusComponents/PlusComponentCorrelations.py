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
        Универсальный метод расчета любого свойства
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

        self.data['Kw'] = self.calculate_property('Kw')

        self.data['critical_temperature'] = self.calculate_property('critical_temperature')

        self.data['critical_pressure'] = self.calculate_property('critical_pressure')

        self.data['acentric_factor'] = self.calculate_property('acentric_factor')

        self.data['critical_volume'] = self.calculate_property('critical_volume')

        self.data['shift_parameter'] = self.calculate_property('shift_parameter')
