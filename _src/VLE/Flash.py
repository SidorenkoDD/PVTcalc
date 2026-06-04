from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton
from _src.Utils.Conditions import Conditions
from _src.Utils.FluidPropertiesCalculatorV2 import FluidPropertiesCalculator
from _src.VLE.FlashResult import FlashResult, PhaseState
from dataclasses import dataclass
from typing import Dict, Any

# (Вставьте сюда классы PhaseState и FlashResult из Шага 1)

class Flash:
    def __init__(self, composition_object: Composition, conditions_object: Conditions):
        self.composition = composition_object
        self.conditions = conditions_object
        self.composition.T = conditions_object.t

    def calculate(self) -> FlashResult: # <-- Указываем тип возврата
        phase_stability_object = TwoPhaseStabilityTest(self.composition, self.conditions.p, self.conditions.t)
        phase_stability_object.calculate_phase_stability()

        if phase_stability_object.stable == False:
            # === ДВУХФАЗНАЯ СИСТЕМА ===
            phase_equil_object = PhaseEquilibriumNewton(
                self.composition,
                self.conditions.p,
                self.conditions.t,
                phase_stability_object.k_values_for_flash
            )
            self.phase_equil_result = phase_equil_object.find_solve_loop()

            # ВАЖНО: Замените 'V' на реальный ключ мольной доли пара из вашего phase_equil_result
            vapor_frac = self.phase_equil_result['Fv']
            liquid_frac = 1.0 - vapor_frac

            liquid_phase_props_object = FluidPropertiesCalculator(
                self.phase_equil_result['xi_l'], self.composition.composition_data,
                phase_equil_object.eos_liquid, self.conditions.p, self.conditions.t
            )
            self.liquid_phase_props = liquid_phase_props_object.calc_all_properties()

            vapour_phase_props_object = FluidPropertiesCalculator(
                self.phase_equil_result['yi_v'], self.composition.composition_data,
                phase_equil_object.eos_vapour, self.conditions.p, self.conditions.t
            )
            self.vapour_phase_props = vapour_phase_props_object.calc_all_properties()

            return FlashResult(
                vapor=PhaseState(mole_fraction=vapor_frac, composition=self.phase_equil_result['yi_v'], properties=self.vapour_phase_props),
                liquid=PhaseState(mole_fraction=liquid_frac, composition=self.phase_equil_result['xi_l'], properties=self.liquid_phase_props),
                is_two_phase=True
            )

        else:
            # === ОДНОФАЗНАЯ СИСТЕМА (Термодинамический трюк) ===
            one_phase_stability_props_object = FluidPropertiesCalculator(
                self.composition.composition, self.composition.composition_data,
                phase_stability_object.liquid_eos, self.conditions.p, self.conditions.t
            )
            self.one_phase_stability_props = one_phase_stability_props_object.calc_all_properties()

            # ТРЮК: Мы присваиваем весь состав "жидкости" (доля 1.0), 
            # а "пару" даем долю 0.0 и тот же состав. 
            # Для DLE это математически идеально: газа выделилось 0, осталась вся жидкость.
            # (Если нужно определить, газ это или жидкость, можно добавить простую проверку плотности)
            return FlashResult(
                vapor=PhaseState(mole_fraction=0.0, composition=self.composition.composition, properties={}),
                liquid=PhaseState(mole_fraction=1.0, composition=self.composition.composition, properties=self.one_phase_stability_props),
                is_two_phase=False
            )