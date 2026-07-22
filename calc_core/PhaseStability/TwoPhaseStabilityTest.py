"""
Тест термодинамической стабильности Михельсена.

Отвечает на вопрос "распадётся ли фид на две фазы при заданных P/T, или
останется однофазным" — первый шаг любого флэша (см. `VLE.Flash`). Запускает
две независимые пробные итерации (`_loop_vapour`/`_loop_liquid`, стартующие
с K-значений Вильсона) и по их сходимости (`_interpetate_stability_analysis`)
решает: `stable=True` (однофазная система) или `stable=False` +
`k_values_for_flash` (стартовые K для `VLE.PhaseEquilibriumNewton`).
"""

import logging

import numpy as np

from calc_core.Composition.Composition import Composition
from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.Utils.EngineConfig import EngineConfig
from calc_core.Utils.Errors import StopIterationError

logger = logging.getLogger(__name__)


class TwoPhaseStabilityTest:
    """
    Один тест стабильности для заданного состава при заданных P/T.

    После `calculate_phase_stability()`: `self.stable` (bool),
    `self.S_v`/`self.S_l` (суммы пробных мольных чисел — диагностика,
    используется в т.ч. `PhaseIdentificator`), `self.k_values_for_flash`
    (стартовые K для двухфазного расчёта, `None` если стабильно).
    """

    def __init__(
        self, composition: Composition, p: float, t: float,
        *, config: EngineConfig | None = None,
    ):
        """
        Parameters
        ----------
        composition : Composition
            Состав фида (с уже посчитанными `composition_data` —
            `acentric_factor`, `critical_pressure`, `critical_temperature`
            как минимум).
        p : float
            Давление, бар.
        t : float
            Температура, K.
        """

        self.p = float(p)
        self.t = float(t)
        self.config = config or EngineConfig.defaults()

        self.stable = None
        self.k_values_liquid = None
        self.k_values_vapour = None
        self.k_values_for_flash = None

        self.S_l = None
        self.S_v = None

        self.convergence_l = False
        self.convergence_v = False
        self.convergence_trivial_solution_l = False
        self.convergence_trivial_solution_v = False

        self._composition = composition
        self.zi = composition.composition
        self.composition_data = composition.composition_data
        # Векторное внутреннее представление
        self._components = tuple(self.zi.keys())
        self._component_index = {comp: i for i, comp in enumerate(self._components)}
        self._nc = len(self._components)

        self._z_feed = np.fromiter(
            (self.zi[c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )

        self._acentric_factor = np.fromiter(
            (self.composition_data['acentric_factor'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )
        self._critical_pressure = np.fromiter(
            (self.composition_data['critical_pressure'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )
        self._critical_temperature = np.fromiter(
            (self.composition_data['critical_temperature'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )

        self._k_v_arr = None
        self._k_l_arr = None
        self._k_flash_arr = None

        self._yi_v_arr = None
        self._xi_l_arr = None

        self._mixture_fugacities_arr = None

        # Публичные dict-представления для совместимости
        self.yi_v = None
        self.xi_l = None

    # =====================================================================================
    # ВСПОМОГАТЕЛЬНЫЕ МЕТОДЫ
    # =====================================================================================

    def _array_to_dict(self, arr: np.ndarray):
        """Конвертирует numpy-массив (по индексам `self._components`) в `{компонент: значение}`."""
        return {comp: float(arr[i]) for i, comp in enumerate(self._components)}

    def _dict_to_array(self, data: dict):
        """Конвертирует `{компонент: значение}` в numpy-массив (порядок — `self._components`)."""
        return np.fromiter(
            (data[c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )

    @staticmethod
    def _safe_log_metric(arr: np.ndarray):
        """Метрика `sum(ln(arr)^2)` — вспомогательная (используется концептуально как `_check_trivial_solution`, но не вызывается напрямую в текущем коде)."""
        return float(np.sum(np.log(arr) ** 2))

    # =====================================================================================
    # WILSON K-VALUES
    # =====================================================================================

    def _calc_k_values_wilson_array(self):
        """
        Начальное приближение констант равновесия по корреляции Вильсона:
        `K_i = (Pc_i/P) * exp(5.37 (1+ω_i)(1 - Tc_i/T))`. Используется как
        общий старт для обоих пробных циклов (`_loop_vapour`/`_loop_liquid`).

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        return np.exp(
            5.37 * (1.0 + self._acentric_factor) * (1.0 - (self._critical_temperature / self.t))
        ) / (self.p / self._critical_pressure)

    # =====================================================================================
    # РАСЧЕТ ПРОБНЫХ СОСТАВОВ
    # =====================================================================================

    def _calc_Yi_v(self, k_values_vapour: np.ndarray):
        """Ненормированные пробные мольные числа паровой фазы: `Y_i = z_i * K_i`."""
        return self._z_feed * k_values_vapour

    def _calc_Xi_l(self, k_values_liquid: np.ndarray):
        """Ненормированные пробные мольные числа жидкой фазы: `X_i = z_i / K_i`."""
        return self._z_feed / k_values_liquid

    @staticmethod
    def _calc_S(arr: np.ndarray):
        """Сумма пробных мольных чисел (`S_v` или `S_l`) — диагностика стабильности: `S <= 1` указывает на тривиальное/стабильное направление."""
        return float(np.sum(arr))

    @staticmethod
    def _normalize_mole_fractions(arr: np.ndarray, S: float):
        """Нормирует пробные мольные числа на их сумму `S`, получая мольные доли пробной фазы."""
        return arr / S

    # =====================================================================================
    # Ri
    # =====================================================================================

    @staticmethod
    def _calc_ri_vapour(vapour_fugacities: np.ndarray, mixture_fugacities: np.ndarray, S_v: float):
        """Отношение фугитивностей для пробной паровой фазы: `r_i = exp(ln f_mix - ln f_vapour) / S_v`."""
        # ri = exp(ln f_mix - ln f_vapour) / S_v
        return np.exp(mixture_fugacities - vapour_fugacities) / S_v

    @staticmethod
    def _calc_ri_liquid(liquid_fugacities: np.ndarray, mixture_fugacities: np.ndarray, S_l: float):
        """Отношение фугитивностей для пробной жидкой фазы: `r_i = exp(ln f_liquid - ln f_mix) * S_l`."""
        # ri = exp(ln f_liquid - ln f_mix) * S_l
        return np.exp(liquid_fugacities - mixture_fugacities) * S_l

    @staticmethod
    def _update_k_values(k_values_old: np.ndarray, ri: np.ndarray):
        """Итерационное обновление K-значений в методе последовательных подстановок: `K_new = K_old * r_i`."""
        return k_values_old * ri

    def _check_convergence(self, ri: np.ndarray) -> bool:
        """
        Критерий сходимости пробного цикла: `sum((r_i - 1)^2) < TOL_TWO_PHASE_STABILITY_CONVERGENCE`.

        Parameters
        ----------
        ri : np.ndarray

        Returns
        -------
        bool
        """
        sum_sq_ri = float(np.sum((ri - 1.0) ** 2))
        return sum_sq_ri < self.config.stability_convergence_tolerance

    def _check_trivial_solution(self, k_values: np.ndarray):
        """
        Критерий тривиального решения пробного цикла: `sum(ln(K_i)^2) < TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION`,
        то есть пробная фаза сходится к составу самого фида (K_i -> 1).

        Parameters
        ----------
        k_values : np.ndarray

        Returns
        -------
        bool
        """
        return float(np.sum(np.log(k_values) ** 2)) < self.config.stability_trivial_tolerance

    # =====================================================================================
    # ИНТЕРПРЕТАЦИЯ РЕЗУЛЬТАТОВ
    # =====================================================================================

    def _interpetate_stability_analysis(self):
        """
        Классический критерий Михельсена по итогам обоих пробных циклов
        (сходимость к нетривиальному решению + `S_v`/`S_l` относительно 1.0):
        заполняет `self.stable` и, если нестабильно, `self.k_values_for_flash`
        (комбинация `k_values_vapour`/`k_values_liquid` в зависимости от того,
        какая из пробных фаз "сработала").

        Вызывается автоматически в конце `calculate_phase_stability()`.
        """
        self.S_v_rounded = round(self.S_v, 7)
        self.S_l_rounded = round(self.S_l, 7)

        if (
            (self.convergence_trivial_solution_v and self.convergence_trivial_solution_l)
            or ((self.S_v_rounded <= 1.0) and self.convergence_trivial_solution_l)
            or (self.convergence_trivial_solution_v and (self.S_l_rounded <= 1.0))
            or ((self.S_v_rounded <= 1.0) and (self.S_l_rounded <= 1.0))
        ):
            self.stable = True
            self._k_flash_arr = None
            self.k_values_for_flash = None

        elif (self.S_v_rounded > 1.0) and self.convergence_trivial_solution_l:
            self.stable = False
            self._k_flash_arr = self._k_v_arr.copy()
            self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)

        elif self.convergence_trivial_solution_v and (self.S_l_rounded > 1.0):
            self.stable = False
            self._k_flash_arr = self._k_l_arr.copy()
            self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)

        elif (self.S_v_rounded > 1.0) and (self.S_l_rounded > 1.0):
            self.stable = False
            self._k_flash_arr = self._k_v_arr * self._k_l_arr
            self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)

        elif (self.S_v_rounded > 1.0) and (self.S_l_rounded <= 1.0):
            self.stable = False
            self._k_flash_arr = self._k_v_arr.copy()
            self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)

        else:
            self.stable = False
            self._k_flash_arr = self._k_l_arr.copy()
            self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)

    # =====================================================================================
    # ОСНОВНОЙ РАСЧЕТ
    # =====================================================================================

    def calculate_phase_stability(self):
        """
        Главная точка входа: полный тест стабильности Михельсена.

        Считает фугитивности исходного (фидового) состава как точку
        отсчёта, генерирует стартовые K по Вильсону
        (`_calc_k_values_wilson_array`), запускает пробную паровую и
        пробную жидкую итерацию (`_loop_vapour`/`_loop_liquid`) и по их
        результату определяет стабильность (`_interpetate_stability_analysis`).
        После вызова читайте `self.stable`/`self.k_values_for_flash`.
        """
        self.initial_eos = BrusilovskiyEOS(composition=self._composition, p=self.p, t=self.t)
        self.initial_eos.calc_eos()
        self._mixture_fugacities_arr = self.initial_eos.fugacities.copy()

        k_init = self._calc_k_values_wilson_array()
        self._loop_vapour(k_init.copy())
        self._loop_liquid(k_init.copy())

        self._interpetate_stability_analysis()
        logger.debug(
            "P=%s, T=%s: stable=%s, S_v=%s, S_l=%s",
            self.p, self.t, self.stable, self.S_v, self.S_l,
        )

    # =====================================================================================
    # ЦИКЛ ПО VAPOUR-ТЕСТУ
    # =====================================================================================

    def _loop_vapour(self, k_values: np.ndarray):
        """
        Пробная итерация "а что если появится паровая фаза" — метод
        последовательных подстановок по фугитивностям (не Ньютон, в отличие
        от `VLE.PhaseEquilibriumNewton`). Заполняет `self.S_v`,
        `self.k_values_vapour`, `self.convergence_v`/`self.convergence_trivial_solution_v`.

        Parameters
        ----------
        k_values : np.ndarray
            Стартовое приближение K (обычно из `_calc_k_values_wilson_array`).

        Raises
        ------
        StopIterationError
            Если не сошлось за 100000 итераций.
        """
        i = 0
        self._k_v_arr = k_values.copy()
        self.k_values_vapour = self._array_to_dict(self._k_v_arr)

        while True:
            Yi_v = self._calc_Yi_v(self._k_v_arr)
            self.S_v = self._calc_S(Yi_v)
            self._yi_v_arr = self._normalize_mole_fractions(Yi_v, self.S_v)
            self.yi_v = self._array_to_dict(self._yi_v_arr)

            new_composition = self._composition.new_composition(self.yi_v)
            self.vapour_eos = BrusilovskiyEOS(composition=new_composition, p=self.p, t=self.t)
            self.vapour_eos.calc_eos()
            vapour_fugacities = self.vapour_eos.fugacities

            ri_v = self._calc_ri_vapour(
                vapour_fugacities=vapour_fugacities,
                mixture_fugacities=self._mixture_fugacities_arr,
                S_v=self.S_v,
            )

            if self._check_convergence(ri_v):
                self.convergence_v = True
                break

            self._k_v_arr = self._update_k_values(self._k_v_arr, ri_v)
            self.k_values_vapour = self._array_to_dict(self._k_v_arr)

            if self._check_trivial_solution(self._k_v_arr):
                self.convergence_trivial_solution_v = True
                break

            i += 1
            if i > self.config.stability_max_iterations:
                logger.warning(
                    "_loop_vapour: не сошёлся за %s итераций (P=%s, T=%s)",
                    self.config.stability_max_iterations, self.p, self.t,
                )
                raise StopIterationError(
                    f'Число итераций теста стабильности превысило '
                    f'{self.config.stability_max_iterations}'
                )

    # =====================================================================================
    # ЦИКЛ ПО LIQUID-ТЕСТУ
    # =====================================================================================

    def _loop_liquid(self, k_values: np.ndarray):
        """
        Пробная итерация "а что если появится жидкая фаза" — зеркало
        `_loop_vapour` для жидкости. Заполняет `self.S_l`,
        `self.k_values_liquid`, `self.convergence_l`/`self.convergence_trivial_solution_l`.

        Parameters
        ----------
        k_values : np.ndarray
            Стартовое приближение K (обычно из `_calc_k_values_wilson_array`).

        Raises
        ------
        StopIterationError
            Если не сошлось за 100000 итераций.
        """
        i = 0
        self._k_l_arr = k_values.copy()
        self.k_values_liquid = self._array_to_dict(self._k_l_arr)

        while True:
            Xi_l = self._calc_Xi_l(self._k_l_arr)
            self.S_l = self._calc_S(Xi_l)
            self._xi_l_arr = self._normalize_mole_fractions(Xi_l, self.S_l)
            self.xi_l = self._array_to_dict(self._xi_l_arr)

            new_composition = self._composition.new_composition(self.xi_l)
            self.liquid_eos = BrusilovskiyEOS(composition=new_composition, p=self.p, t=self.t)
            self.liquid_eos.calc_eos()
            liquid_fugacities = self.liquid_eos.fugacities

            ri_l = self._calc_ri_liquid(
                liquid_fugacities=liquid_fugacities,
                mixture_fugacities=self._mixture_fugacities_arr,
                S_l=self.S_l,
            )

            if self._check_convergence(ri_l):
                self.convergence_l = True
                break

            self._k_l_arr = self._update_k_values(self._k_l_arr, ri_l)
            self.k_values_liquid = self._array_to_dict(self._k_l_arr)

            if self._check_trivial_solution(self._k_l_arr):
                self.convergence_trivial_solution_l = True
                break

            i += 1
            if i > self.config.stability_max_iterations:
                logger.warning(
                    "_loop_liquid: не сошёлся за %s итераций (P=%s, T=%s)",
                    self.config.stability_max_iterations, self.p, self.t,
                )
                raise StopIterationError(
                    f'Число итераций теста стабильности превысило '
                    f'{self.config.stability_max_iterations}'
                )

    # =====================================================================================
    # СВОЙСТВО ДЛЯ РАСЧЕТА ТОЧКИ НАСЫЩЕНИЯ
    # =====================================================================================

    @property
    def k_vals_for_sat_point_calculation(self):
        """
        dict | None: K-значения для использования в поиске давления насыщения
        (`PhaseEnvelope/new_methodv2.py::SaturationPressure` и подобные).

        `None`, если тест ещё не запускался или система стабильна. Если
        обе пробные суммы `S_v`/`S_l` > 1 — берёт K от той пробной фазы,
        у которой сумма больше (эвристика "более выраженное" направление
        нестабильности); иначе — то же, что `k_values_for_flash`.
        """
        if (self.stable is None) or self.stable:
            return None

        if (self.S_v > 1.0) and (self.S_l > 1.0):
            if self.S_v > self.S_l:
                return self.k_values_vapour
            return self.k_values_liquid

        return self.k_values_for_flash
