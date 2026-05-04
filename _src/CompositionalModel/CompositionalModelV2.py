from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton 
class CompositionalModel:

    def __init__(self, composition : Composition):
        self.composition = composition

    def flash(self, P, T):
        phase_stability_obj = TwoPhaseStabilityTest(composition = self.composition,
                                                    p = P, 
                                                    t = T)
        phase_stability_obj.calculate_phase_stability()
        phase_equil_obj = PhaseEquilibriumNewton(composition=self.composition,
                                                 p = P,
                                                 t = T,
                                                 k_values = phase_stability_obj.k_values_for_flash)
        
        result = phase_equil_obj.find_solve_loop()
        return result