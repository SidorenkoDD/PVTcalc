"""Единый публичный фасад расчётов над составом флюида.

CompositionalModel связывает состав с фасадами Flash, экспериментов и
фазовой огибающей. Все публичные методы принимают инженерные единицы,
документированные в сигнатурах, и не мутируют исходный Composition.
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.Experiments.ExperimentsFacade import ExperimentsFacade
from calc_core.PhaseEnvelope.PhaseEnvelopeFacade import PhaseEnvelopeFacade
from calc_core.Utils.Conditions import Conditions
from calc_core.Utils.Cancellation import CancellationToken, ProgressCallback, report_progress
from calc_core.Utils.EngineConfig import EngineConfig
from calc_core.Utils.Results import ResultStore
from calc_core.Utils.Validation import (
    validate_composition_normalized,
    validate_positive_pressure,
    validate_temperature_celsius,
)
from calc_core.VLE.Flash import Flash
from calc_core.VLE.FlashResult import FlashResult

logger = logging.getLogger(__name__)


class CompositionalModel:
    """Обёртка над составом с историей расчётов."""

    def __init__(self, composition: Composition, config: EngineConfig | None = None):
        """
        Parameters
        ----------
        composition : Composition
            Состав с уже посчитанными composition_data и нормированными
            мольными долями (сумма = 1).

        Raises
        ------
        InputValidationError
            Если сумма мольных долей состава != 1.
        """
        validate_composition_normalized(composition.composition)
        self.composition = composition
        self.config = config or EngineConfig.defaults()
        self.result_store_object = ResultStore()
        self.phase_envelope = PhaseEnvelopeFacade(self)
        self.experiments = ExperimentsFacade(self)

    def flash(
        self,
        p_bar: float,
        t_celsius: float,
        *,
        cancellation_token: CancellationToken | None = None,
        progress_callback: ProgressCallback | None = None,
    ) -> FlashResult:
        """Рассчитывает одно- или двухфазный Flash через единый контракт.

        Parameters
        ----------
        p_bar : float
            Давление, bar.
        t_celsius : float
            Температура, °C. Внутри расчётного ядра она переводится в Kelvin.

        Returns
        -------
        FlashResult
            Типизированный результат с обеими фазами, phase_type для
            однофазного состояния и физическими diagnostics.

        Notes
        -----
        Flash меняет температуру и EOS-зависимые свойства переданного
        Composition, поэтому расчёт выполняется на независимой глубокой
        копии. Исходный состав модели остаётся неизменным.
        """
        validate_positive_pressure(p_bar, name="p_bar")
        validate_temperature_celsius(t_celsius, name="t_celsius")
        if cancellation_token is not None:
            cancellation_token.throw_if_cancelled()
        report_progress(progress_callback, 0.0, "Preparing flash")
        logger.info(
            "CompositionalModel.flash: P=%s bar, T=%s °C", p_bar, t_celsius,
        )

        work = self.composition.new_composition(
            self.composition.composition,
            deep_copy=True,
        )
        result = Flash(
            work, Conditions(p_bar, t_celsius), config=self.config,
            cancellation_token=cancellation_token,
            progress_callback=progress_callback,
        ).calculate()
        report_progress(progress_callback, 1.0, "Flash complete")
        self.result_store_object.add(
            module="Flash",
            params={"p_bar": p_bar, "t_celsius": t_celsius},
            data=result,
        )
        return result
