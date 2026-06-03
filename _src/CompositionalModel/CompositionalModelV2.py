from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton
from _src.Utils.Results2 import ResultStore


class CompositionalModel:

    def __init__(self, composition : Composition):
        self.composition = composition
        self.result_store_object = ResultStore()

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

        self.result_store_object.add(module = 'Flash', params= {'P': P, 'T':T}, data= result)
        return result