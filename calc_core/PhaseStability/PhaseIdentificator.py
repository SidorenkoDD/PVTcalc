from __future__ import annotations

import logging

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from calc_core.VLE.Flash import Flash

logger = logging.getLogger(__name__)


class PhaseIdentificator:
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
            logger.debug('Sl: %s', (self.flash_object.phase_stability_object.S_l, self.flash_object.phase_stability_object.S_l_rounded))
            logger.debug('Sv: %s', (self.flash_object.phase_stability_object.S_v, self.flash_object.phase_stability_object.S_v_rounded))
            logger.debug('Sv > Sl : %s', self.flash_object.phase_stability_object.S_v > self.flash_object.phase_stability_object.S_l)
            logger.debug('Sv > Sl rounded: %s', self.flash_object.phase_stability_object.S_v_rounded > self.flash_object.phase_stability_object.S_l_rounded)

            stability_phase = 'vapor' if Sv < Sl else 'liquid'
            volume_phase = 'vapor' if res > self.Constant * 10 else 'liquid'
            logger.debug('stability heuristic -> %s', stability_phase)
            logger.debug('volume heuristic (%s) -> %s', res, volume_phase)
            return stability_phase if stability_phase == volume_phase else 'ambiguous'
        return None


# Обратная совместимость с исторической опечаткой.
PhaseIndentificator = PhaseIdentificator
