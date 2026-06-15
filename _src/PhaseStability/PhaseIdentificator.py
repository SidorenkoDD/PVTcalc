from _src.VLE.Flash import Flash
import numpy as np

class PhaseIndentificator:
    def __init__(self, flash_object: Flash, Constant = 1.75):
        self.flash_object = flash_object
        self.Constant = Constant

    def identify_phase(self):
        if self.flash_object.phase_stability_object.stable:
            sum_b = float(self.flash_object.phase_stability_object.liquid_eos._b @ self.flash_object.phase_stability_object.liquid_eos.composition_vector)
            molar_vol = self.flash_object.one_phase_stability_props['molar_volume']

            print(f'b_arr : {self.flash_object.phase_stability_object.liquid_eos._b}')
            print(f'Composition : {self.flash_object.phase_stability_object.liquid_eos.composition_vector}')
            print(f'Sum_b: {sum_b}')
            print(f'Mol_vol: {molar_vol}')
            print(f'Res: {molar_vol / sum_b}')