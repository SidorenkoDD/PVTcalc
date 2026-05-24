from _src.Composition.CompositionV2 import Composition
from _src.Utils.Conditions import Conditions
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton
from _src.Utils.FluidPropertiesCalculatorV2 import FluidPropertiesCalculator


class TwoPhaseFlash:
    def __init__(self, composition : Composition, conditions : Conditions):
        self._composition = composition
        self._conditions = conditions
        self._stability = None

    def _calculate_stability_test(self):
        self._stability_test_object = TwoPhaseStabilityTest(self._composition,
                                                      p = self._conditions.p,
                                                      t = self._conditions.t)
        self._stability_test_object.calculate_phase_stability()
        self._stability = self._stability_test_object.stable
        print(self._stability)

    def calculate_flash(self):
        self._calculate_stability_test()

        if self._stability == False:
            print('TwoPhase flash')
            self._phase_equilibrium_object = PhaseEquilibriumNewton(composition=self._composition,
                                                                   p = self._conditions.p,
                                                                   t = self._conditions.t,
                                                                   k_values = self._stability_test_object.k_values_for_flash)
            self._phase_equilibrium_object.find_solve_loop()
            self._vapour_phase_props = FluidPropertiesCalculator(composition = self._phase_equilibrium_object.yi_v,
                                                                 composition_properties = self._composition._composition_data,
                                                                 eos_object = self._phase_equilibrium_object.eos_vapour,
                                                                 p = self._conditions.p,
                                                                 T = self._conditions.t)
            print(self._vapour_phase_props.density)
            self._liquid_phase_props = FluidPropertiesCalculator(composition = self._phase_equilibrium_object.xi_l,
                                                                 composition_properties = self._composition._composition_data,
                                                                 eos_object = self._phase_equilibrium_object.eos_liquid,
                                                                 p = self._conditions.p,
                                                                 T = self._conditions.t)
            print(self._liquid_phase_props.density)
        else:
            print('one_phase_props')
            self._one_phase_props = FluidPropertiesCalculator(composition = self._phase_equilibrium_object.yi_v,
                                                                 composition_properties = self._composition._composition_data,
                                                                 eos_object = self._phase_equilibrium_object.eos_vapour,
                                                                 p = self._conditions.p,
                                                                 T = self._conditions.t)
            
            
        

