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
            res = molar_vol / sum_b
            Sv = self.flash_object.phase_stability_object.S_v
            Sl = self.flash_object.phase_stability_object.S_l
            print(f'Sl: {self.flash_object.phase_stability_object.S_l, self.flash_object.phase_stability_object.S_l_rounded}')
            print(f'Sv: {self.flash_object.phase_stability_object.S_v, self.flash_object.phase_stability_object.S_v_rounded}')
            print(f'Sv > Sl : {self.flash_object.phase_stability_object.S_v > self.flash_object.phase_stability_object.S_l}')
            print(f'Sv > Sl rounded: {self.flash_object.phase_stability_object.S_v_rounded > self.flash_object.phase_stability_object.S_l_rounded}')

            if Sv < Sl:
                print(f'Sv < Sl -> vapour')
            else:
                print(f'Sv > Sl -> liquid')
            if res > self.Constant * 10:
                print(res)
                print('vapour')
            else:
                print(res)
                print('liquid')