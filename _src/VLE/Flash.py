from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton
from _src.Utils.Conditions import Conditions
from _src.Utils.FluidPropertiesCalculatorV2 import FluidPropertiesCalculator


class Flash:
    def __init__(self, composition_object : Composition, conditions_object : Conditions):
        self.composition = composition_object
        self.conditions = conditions_object
        self.composition.T = conditions_object.t


    def calculate(self):
        phase_stability_object = TwoPhaseStabilityTest(self.composition, self.conditions.p, self.conditions.t)
        phase_stability_object.calculate_phase_stability()


        if phase_stability_object.stable == False:

            phase_equil_object = PhaseEquilibriumNewton(self.composition,
                                                        self.conditions.p,
                                                        self.conditions.t,
                                                        phase_stability_object.k_values_for_flash)

            self.phase_equil_result = phase_equil_object.find_solve_loop()

            liquid_phase_props_object = FluidPropertiesCalculator(self.phase_equil_result['xi_l'],
                                                        self.composition.composition_data,
                                                        phase_equil_object.eos_liquid,
                                                        self.conditions.p,
                                                        self.conditions.t)
            self.liquid_phase_props = liquid_phase_props_object.calc_all_properties()

            vapour_phase_props_object = FluidPropertiesCalculator(self.phase_equil_result['yi_v'],
                                                        self.composition.composition_data,
                                                        phase_equil_object.eos_vapour,
                                                        self.conditions.p,
                                                        self.conditions.t)
            self.vapour_phase_props = vapour_phase_props_object.calc_all_properties()

            return {'phase_equil' : self.phase_equil_result,
                    'liquid_props' : self.liquid_phase_props,
                    'vapour_props' : self.vapour_phase_props}

        else:

            one_phase_stability_props_object = FluidPropertiesCalculator(self.composition.composition,
                                                                       self.composition.composition_data,
                                                                       phase_stability_object.vapour_eos,
                                                                       self.conditions.p,
                                                                       self.conditions.t)
            
            self.one_phase_stability_props = one_phase_stability_props_object.calc_all_properties()

            return {'one_phase_props' : self.one_phase_stability_props}