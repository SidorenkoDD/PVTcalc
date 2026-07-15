"""
Тонкий фасад над стабильность+равновесие с логированием результатов.

Альтернатива `VLE.Flash` для тех же двух шагов (тест стабильности → решение
равновесия), но: (1) не обрабатывает однофазный случай — `flash()` всегда
идёт по пути `PhaseEquilibriumNewton`, даже если система стабильна (см.
caveat в docstring `flash`); (2) вместо `FlashResult` возвращает сырой dict
из `PhaseEquilibriumNewton.find_solve_loop()` и попутно логирует каждый
вызов в `ResultStore`. Используется в `test_notebook.ipynb`
("Работа через интерфейс Model").
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.VLE.PhaseEquilibriumNewton import PhaseEquilibriumNewton
from calc_core.Utils.Results import ResultStore
from calc_core.PhaseEnvelope.PhaseEnvelopeFacade import PhaseEnvelopeFacade

logger = logging.getLogger(__name__)


class CompositionalModel:
    """Обёртка над составом с историей расчётов (`result_store_object`, см. `Utils.Results.ResultStore`)."""

    def __init__(self, composition : Composition):
        """
        Parameters
        ----------
        composition : Composition
            Состав с уже посчитанными `composition_data`.
        """
        self.composition = composition
        self.result_store_object = ResultStore()
        self.phase_envelope = PhaseEnvelopeFacade(self)

    def flash(self, P, T):
        """
        Тест стабильности + решение двухфазного равновесия для заданных P/T,
        результат логируется в `self.result_store_object`.

        **Caveat**: в отличие от `VLE.Flash.calculate()`, здесь нет ветки
        для стабильного (однофазного) случая — если
        `phase_stability_obj.k_values_for_flash` окажется `None` (система
        стабильна), `PhaseEquilibriumNewton` получит на вход `None` вместо
        K-значений и, скорее всего, упадёт. Использовать для точек, заведомо
        находящихся в двухфазной области.

        Parameters
        ----------
        P : float
            Давление, бар.
        T : float
            Температура, K.

        Returns
        -------
        dict
            То же, что возвращает `PhaseEquilibriumNewton.find_solve_loop()`:
            `{"yi_v", "xi_l", "Ki", "Fv", "Fl", "Z_v", "Z_l"}`.
        """
        logger.debug("CompositionalModel.flash: P=%s, T=%s", P, T)
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