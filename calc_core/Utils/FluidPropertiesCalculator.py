"""
Расчёт физических свойств одной фазы по уже посчитанному EOS.

Принимает состав фазы (мольные доли + свойства компонент) и объект
`BrusilovskiyEOS` с уже выбранным Z-фактором, отдаёт молярную массу, объём,
плотность, Z со сдвигом и вязкость (по умолчанию LBC). Используется
`VLE.Flash` (по разу на каждую фазу) и всеми PVT-экспериментами.
"""

from calc_core.Utils.Constants import CONSTANT_R
from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.Utils.Viscosity import ViscosityFactory


class FluidPropertiesCalculator:
    """Ленивый (через `@property`) расчёт свойств одной фазы; ничего не считается до первого обращения к конкретному свойству."""

    def __init__(self, composition: dict, composition_properties: dict, eos_object: BrusilovskiyEOS, p: float, T: float,
                 viscosity_method: str = 'LBC'):
        """
        Parameters
        ----------
        composition : dict
            `{компонент: мольная_доля}` **этой фазы** (не всего фида).
        composition_properties : dict
            `composition_data` фида (`Composition.composition_data`) —
            статические свойства компонент (молярная масса, Kw и т.д.),
            общие для всех фаз.
        eos_object : BrusilovskiyEOS
            Уже посчитанный (`calc_eos()` вызван) объект EOS **для этой
            фазы** (построен на составе именно этой фазы).
        p : float
            Давление, бар.
        T : float
            Температура, K.
        viscosity_method : str, optional
            Ключ метода вязкости для `Viscosity.ViscosityFactory`. По
            умолчанию `'LBC'` (единственный реально реализованный метод —
            `'PEDERSEN'` мапится на тот же класс, см. `Utils/Viscosity.py`).
        """
        self._composition = composition
        self._composition_properties = composition_properties
        self._eos = eos_object
        self._p = p
        self._T = T
        self._viscosity_method = viscosity_method

    @property
    def molar_mass(self):
        """float: молярная масса фазы — средневзвешенное по мольным долям."""
        M = 0.0

        for component in self._composition:
            M += self._composition[component] * self._composition_properties['molar_mass'][component]

        return M

    @property
    def molar_volume(self):
        """float: молярный объём фазы, RT·Z/P со сдвигом Пенелу (`self._eos.shift_parametr`)."""
        return CONSTANT_R * 10 * self._T * self._eos.z / self._p - self._eos.shift_parametr

    @property
    def molar_density(self):
        """float: мольная плотность фазы (`1 / molar_volume`)."""
        return 1 / self.molar_volume

    @property
    def density(self):
        """float: массовая плотность фазы (`molar_density * molar_mass`)."""
        return self.molar_density * self.molar_mass

    @property
    def z_shift(self):
        """float: Z-фактор с поправкой на volume-shift (в отличие от `self._eos.z` — "сырого" Z без сдвига)."""
        return self._eos.z - self._p * self._eos.shift_parametr / (self._T * CONSTANT_R * 10)

    @property
    def viscosity(self):
        """
        float: вязкость фазы, сП — по методу `self._viscosity_method`
        (см. `Utils.Viscosity.ViscosityFactory`). Коэффициенты `alpha0..alpha4`
        сейчас захардкожены (см. закомментированные альтернативные наборы
        ниже в коде — судя по всему, разные попытки калибровки).
        """
        visc = ViscosityFactory().create_viscosity_object(self._viscosity_method)
        return visc(mole_fractions=self._composition, composition_data=self._composition_properties,
                   phase_density=self.density, mw=self.molar_mass, temperature=self._T,
                    alpha0=0.1023,
                    alpha1=0.023364,
                    alpha2=0.058533,
                    alpha3=-0.040758,
                    alpha4=0.0093324,
                    ).calculate()

        # return visc(mole_fractions=self._composition, composition_data=self._composition_properties,
        #             phase_density=self.density, mw=self.molar_mass, temperature=self._T,
        #             alpha0=0.1314433068037033,
        #             alpha1=0.015912599861621857,
        #             alpha2=0.05853300169110298,
        #             alpha3=-0.04075799882411957,
        #             alpha4=0.009332399815320969,
        #             ).calculate()
        # return visc(mole_fractions=self._composition, composition_data=self._composition_properties,
        #             phase_density=self.density, mw=self.molar_mass, temperature=self._T).calculate()
    
    def calc_all_properties(self):
        """
        Считает все свойства фазы разом.

        Returns
        -------
        dict
            Ключи: `'molecular_ weight'` (с пробелом внутри — опечатка в
            коде, см. CLAUDE.md Known Issues), `'molar_volume'`,
            `'molar_density'`, `'density'`, `'z'` (со сдвигом,
            т.е. `z_shift`), `'viscosity'`. Именно этот словарь оказывается
            в `PhaseState.properties` в `FlashResult`.
        """

        return {'molecular_ weight': self.molar_mass,
                'molar_volume': self.molar_volume,
                'molar_density': self.molar_density,
                'density': self.density,
                'z' : self.z_shift,
                'viscosity': self.viscosity}