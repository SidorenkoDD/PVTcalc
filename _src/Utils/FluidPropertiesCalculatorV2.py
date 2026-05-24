from _src.Utils.Constants import CONSTANT_R
from _src.EOS.BrusilovskiyEOSV2 import BrusilovskiyEOS
from _src.Utils.Viscosity import ViscosityFactory


class FluidPropertiesCalculator:
    '''
    Модуль для расчета свойств фазы после флеша

    Args:
        composition: dict с мольными долями состава фазы
        composition_properties: dict со свойствами компонент фазы
        eos_object: объект EOS для доступа к z и шифт параметру
        P: давление в БАР
        t: температура в К
        viscosity_method: метод для расчета вязкости, по умолчанию LBС

    TODO:
        Надо разобраться с тем, почему где то используется 8.31, а где то 83.1
    '''
    def __init__(self, composition: dict, composition_properties: dict, eos_object: BrusilovskiyEOS, p: float, T: float,
                 viscosity_method: str = 'LBC'):
        self._composition = composition
        self._composition_properties = composition_properties
        self._eos = eos_object
        self._p = p
        self._T = T
        self._viscosity_method = viscosity_method

    @property
    def molar_mass(self):
        M = 0.0

        for component in self._composition:
            M += self._composition[component] * self._composition_properties['molar_mass'][component]

        return M

    @property
    def molar_volume(self):
        return (CONSTANT_R * 10 * self._T * self._eos.z / self._p - self._eos.shift_parametr)

    @property
    def molar_density(self):
        return 1 / self.molar_volume

    @property
    def density(self):
        return self.molar_density * self.molar_mass

    @property
    def z_shift(self):
        return self._eos.z - self._p * self._eos.shift_parametr / (self._T * CONSTANT_R * 10)

    @property
    def viscosity(self):
        visc = ViscosityFactory().create_viscosity_object(self._viscosity_method)
        return visc(mole_fractions=self._composition, composition_data=self._composition_properties,
                   phase_density=self.density, mw=self.molar_mass, temperature=self._T,
                    alpha0=0.12538570165634155,
                    alpha1=0.02722247689962387,
                    alpha2=0.004341159015893936,
                    alpha3=-0.004339530132710934,
                    alpha4=0.0007456476450897753,
                    ).calculate()

