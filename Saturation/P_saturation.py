from _src.Composition.Composition import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest


class SaturationPressureCalc:
    def __init__(self, composition : Composition, p_max, temperature, p_min = 0.01):
        self.composition = composition
        self.p_max_init = p_max
        self.temperature = temperature
        self.p_min = p_min

    