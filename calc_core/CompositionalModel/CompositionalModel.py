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
from calc_core.Utils.Validation import validate_composition_normalized, validate_positive_pressure, validate_temperature_kelvin
from calc_core.PhaseEnvelope.PhaseEnvelopeFacade import PhaseEnvelopeFacade
from calc_core.Experiments.ExperimentsFacade import ExperimentsFacade

logger = logging.getLogger(__name__)


class CompositionalModel:
    """Обёртка над составом с историей расчётов (`result_store_object`, см. `Utils.Results.ResultStore`)."""

    def __init__(self, composition : Composition):
        """
        Parameters
        ----------
        composition : Composition
            Состав с уже посчитанными `composition_data` и нормированными
            мольными долями (сумма = 1). Ненормированный состав отклоняется —
            вызовите `composition.normalize_composition()` заранее.

        Raises
        ------
        InputValidationError
            Если сумма мольных долей состава != 1.
        """
        validate_composition_normalized(composition.composition)
        self.composition = composition
        self.result_store_object = ResultStore()
        self.phase_envelope = PhaseEnvelopeFacade(self)
        self.experiments = ExperimentsFacade(self)

    def flash(self, P, T_k):
        """
        Тест стабильности + решение двухфазного равновесия для заданных P/T,
        результат логируется в `self.result_store_object`.

        .. deprecated::
            **Для нового кода используйте `calc_core.VLE.Flash.Flash`.**
            Этот метод — единственный из трёх фасадов (`flash`,
            `phase_envelope`, `experiments`), который не приведён к общему
            контракту, и отличается от них тремя вещами:

            1. **Не копирует состав.** `PhaseEquilibriumNewton` мутирует
               `self.composition` через тяжёлый setter `Composition.T`, так
               что `model.composition` тихо меняется после каждого вызова.
               `PhaseEnvelopeFacade`/`ExperimentsFacade` передают в калькулятор
               независимую копию (`_composition_copy`), этот метод — нет.
            2. **Нет ветки для однофазного случая.** Если система стабильна,
               `phase_stability_obj.k_values_for_flash` окажется `None`,
               `PhaseEquilibriumNewton` получит `None` вместо K-значений и,
               скорее всего, упадёт. Использовать только для точек, заведомо
               находящихся в двухфазной области.
            3. **Возвращает сырой dict, а не `FlashResult`** — без свойств фаз
               (плотность, Z, вязкость), которые считает `VLE.Flash`.

            GUI по этим причинам вызывает `VLE.Flash` напрямую
            (см. `gui/services/flash_service.py`). Метод оставлен рабочим ради
            `test_notebook.ipynb` и обратной совместимости.

        Parameters
        ----------
        P : float
            Давление, бар.
        T_k : float
            Температура, **K**.

        Returns
        -------
        dict
            То же, что возвращает `PhaseEquilibriumNewton.find_solve_loop()`:
            `{"yi_v", "xi_l", "Ki", "Fv", "Fl", "Z_v", "Z_l"}`.

        Raises
        ------
        InputValidationError
            Если `P` <= 0 или `T_k` <= 0 K.
        """
        validate_positive_pressure(P, name='P')
        validate_temperature_kelvin(T_k, name='T_k')
        logger.debug("CompositionalModel.flash: P=%s, T=%s K", P, T_k)
        phase_stability_obj = TwoPhaseStabilityTest(composition = self.composition,
                                                    p = P,
                                                    t = T_k)
        phase_stability_obj.calculate_phase_stability()
        phase_equil_obj = PhaseEquilibriumNewton(composition=self.composition,
                                                 p = P,
                                                 t = T_k,
                                                 k_values = phase_stability_obj.k_values_for_flash)

        result = phase_equil_obj.find_solve_loop()

        self.result_store_object.add(module = 'Flash', params= {'P': P, 'T_k': T_k}, data= result)
        return result