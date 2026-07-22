"""
Оркестратор одного двухфазного флэш-расчёта.

Верхнеуровневый узел расчётного пайплайна: принимает состав (`Composition`)
и термобарические условия (`Conditions`), сначала прогоняет тест
термодинамической стабильности (`TwoPhaseStabilityTest`), затем — в
зависимости от результата — либо решает двухфазное равновесие
(`PhaseEquilibriumNewton`), либо считает свойства всего состава как
однофазной жидкости ("ТРЮК", см. docstring `Flash.calculate`). Используется
напрямую в `test_notebook.ipynb` и как строительный блок в
`Experiments/DLE.py`, `CCE.py`, `SeparatorTest.py` (там `Flash` вызывается
многократно, по одному разу на ступень) и в `CompositionalModel`.
"""

import logging
from dataclasses import replace

from calc_core.Composition.Composition import Composition
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.Utils.Conditions import Conditions
from calc_core.Utils.Cancellation import (
    CancellationToken,
    ProgressCallback,
    report_progress,
)
from calc_core.Utils.EngineConfig import EngineConfig
from calc_core.Utils.FluidPropertiesCalculator import FluidPropertiesCalculator
from calc_core.Utils.ResultDiagnostics import ResultWarning, diagnose_flash_result
from calc_core.VLE.FlashResult import FlashResult, PhaseState
from calc_core.VLE.PhaseEquilibriumNewton import PhaseEquilibriumNewton

logger = logging.getLogger(__name__)

class Flash:
    """
    Один флэш-расчёт для заданного состава при заданных P/T.

    Экземпляр — одноразовый: создаётся под конкретную пару (`Composition`,
    `Conditions`), после `calculate()` держит промежуточные объекты расчёта
    как атрибуты (`phase_stability_object`, `phase_equil_result` и т.п.) —
    удобно для отладки, но не предполагает повторного вызова `calculate()`
    с другими условиями (для этого создавайте новый `Flash`).
    """

    def __init__(
        self,
        composition_object: Composition,
        conditions_object: Conditions,
        *,
        config: EngineConfig | None = None,
        cancellation_token: CancellationToken | None = None,
        progress_callback: ProgressCallback | None = None,
    ):
        """
        Parameters
        ----------
        composition_object : Composition
            Состав флюида. Побочный эффект: `composition_object.T` сразу
            перезаписывается температурой из `conditions_object` (Кельвины).
        conditions_object : Conditions
            Давление (бар) и температура (K, уже сконвертированная из °C).
        """
        self.composition = composition_object
        self.conditions = conditions_object
        self.config = config or EngineConfig.defaults()
        self.cancellation_token = cancellation_token
        self.progress_callback = progress_callback
        self.composition.T = conditions_object.t

    def calculate(self) -> FlashResult: # <-- Указываем тип возврата
        """
        Выполняет тест стабильности и в зависимости от результата — двухфазный
        расчёт равновесия либо однофазный "трюк".

        Ветка **двухфазная** (`phase_stability_object.stable is False`):
        решает равновесие Ньютоном по фугитивностям (`PhaseEquilibriumNewton`),
        считает свойства паровой и жидкой фазы отдельно
        (`FluidPropertiesCalculator`).

        Ветка **однофазная** ("ТРЮК"): без проверки, газ это или жидкость,
        весь состав жёстко приписывается жидкости (`liquid.mole_fraction = 1.0`,
        `vapor.mole_fraction = 0.0`, `vapor.properties = {}`). Корректно для
        DLE (газа не выделилось), но не отвечает на вопрос "что это за фаза" —
        для этого предназначен (пока не интегрирован сюда)
        `PhaseStability/PhaseIdentificator.py`.

        Returns
        -------
        FlashResult
            `pressure`, `temperature`, `vapor: PhaseState`, `liquid: PhaseState`,
            `is_two_phase: bool`. `PhaseState.properties` — словарь с ключами
            `molecular_ weight` (со пробелом — опечатка в коде, см. CLAUDE.md),
            `molar_volume`, `molar_density`, `density`, `z`, `viscosity`.
        """
        logger.debug("Flash.calculate(): P=%s бар, T=%s K", self.conditions.p, self.conditions.t)
        if self.cancellation_token is not None:
            self.cancellation_token.throw_if_cancelled()
        report_progress(self.progress_callback, 0.05, "Testing phase stability")

        self.phase_stability_object = TwoPhaseStabilityTest(
            self.composition, self.conditions.p, self.conditions.t,
            config=self.config,
            cancellation_token=self.cancellation_token,
            progress_callback=self.progress_callback,
        )
        self.phase_stability_object.calculate_phase_stability()

        if not self.phase_stability_object.stable:
            logger.debug("P=%s, T=%s: система нестабильна, двухфазный расчёт", self.conditions.p, self.conditions.t)
            # === ДВУХФАЗНАЯ СИСТЕМА ===
            phase_equil_object = PhaseEquilibriumNewton(
                self.composition,
                self.conditions.p,
                self.conditions.t,
                self.phase_stability_object.k_values_for_flash,
                config=self.config,
                cancellation_token=self.cancellation_token,
                progress_callback=self.progress_callback,
            )
            self.phase_equil_result = phase_equil_object.find_solve_loop()
            report_progress(self.progress_callback, 0.82, "Calculating phase properties")

            # ВАЖНО: Замените 'V' на реальный ключ мольной доли пара из вашего phase_equil_result
            vapor_frac = self.phase_equil_result['Fv']
            liquid_frac = 1.0 - vapor_frac
            logger.debug("Flash: сошлось, Fv=%.6f", vapor_frac)

            liquid_phase_props_object = FluidPropertiesCalculator(
                self.phase_equil_result['xi_l'], self.composition.composition_data,
                phase_equil_object.eos_liquid, self.conditions.p, self.conditions.t
            )
            self.liquid_phase_props = liquid_phase_props_object.calc_all_properties()

            vapour_phase_props_object = FluidPropertiesCalculator(
                self.phase_equil_result['yi_v'], self.composition.composition_data,
                phase_equil_object.eos_vapour, self.conditions.p, self.conditions.t
            )
            self.vapour_phase_props = vapour_phase_props_object.calc_all_properties()

            result = FlashResult(pressure = self.conditions.p, temperature = self.conditions.t,
                vapor=PhaseState(mole_fraction=vapor_frac, composition=self.phase_equil_result['yi_v'], properties=self.vapour_phase_props),
                liquid=PhaseState(mole_fraction=liquid_frac, composition=self.phase_equil_result['xi_l'], properties=self.liquid_phase_props),
                is_two_phase=True,
                engine_config=self.config,
            )
            extra = ()
            if phase_equil_object.trivial_solution:
                extra = (ResultWarning(
                    "trivial_equilibrium",
                    "Итерации пришли к K≈1; двухфазное представление требует проверки.",
                    "equilibrium",
                ),)
            diagnostics = diagnose_flash_result(
                result, self.composition.composition, extra_warnings=extra,
            )
            return replace(result, diagnostics=diagnostics)

        else:
            logger.debug(
                "P=%s, T=%s: система стабильна, весь состав отнесён к жидкости (\"ТРЮК\")",
                self.conditions.p, self.conditions.t,
            )
            # === ОДНОФАЗНАЯ СИСТЕМА (Термодинамический трюк) ===
            self.one_phase_stability_props_object = FluidPropertiesCalculator(
                self.composition.composition, self.composition.composition_data,
                self.phase_stability_object.liquid_eos, self.conditions.p, self.conditions.t
            )
            self.one_phase_stability_props = self.one_phase_stability_props_object.calc_all_properties()
            report_progress(self.progress_callback, 0.82, "Calculating phase properties")
            from calc_core.PhaseStability.PhaseIdentificator import PhaseIdentificator
            phase_type = PhaseIdentificator(self).identify_phase()

            # ТРЮК: Мы присваиваем весь состав "жидкости" (доля 1.0), 
            # а "пару" даем долю 0.0 и тот же состав. 
            # Для DLE это математически идеально: газа выделилось 0, осталась вся жидкость.
            # (Если нужно определить, газ это или жидкость, можно добавить простую проверку плотности)
            result = FlashResult(pressure = self.conditions.p, temperature = self.conditions.t,
                vapor=PhaseState(mole_fraction=0.0, composition=self.composition.composition, properties={}),
                liquid=PhaseState(mole_fraction=1.0, composition=self.composition.composition, properties=self.one_phase_stability_props),
                is_two_phase=False,
                phase_type=phase_type,
                engine_config=self.config,
            )
            return replace(
                result,
                diagnostics=diagnose_flash_result(result, self.composition.composition),
            )
