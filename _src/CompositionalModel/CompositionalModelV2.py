from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton
from _src.Utils.Conditions import Conditions
from _src.VLE.flashV2 import TwoPhaseFlash


class CompositionalModel:

    def __init__(self, composition : Composition):
        self.composition = composition

    def flashV2(self, conditions : Conditions):
        flash_object = TwoPhaseFlash(self.composition, conditions = conditions)
        flash_object.calculate_flash()

    def flash(self, conditions: Conditions):
        phase_stability_obj = TwoPhaseStabilityTest(composition = self.composition,
                                                    p = conditions.p, 
                                                    t = conditions.t)
        phase_stability_obj.calculate_phase_stability()
        phase_equil_obj = PhaseEquilibriumNewton(composition=self.composition,
                                                 p = conditions.p,
                                                 t = conditions.t,
                                                 k_values = phase_stability_obj.k_values_for_flash)
        
        result = phase_equil_obj.find_solve_loop()
        return result
    