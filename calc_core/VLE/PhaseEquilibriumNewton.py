"""
Решатель двухфазного равновесия: Rachford-Rice + Ньютон по фугитивностям.

Принимает состав, условия и начальное приближение K-значений (обычно из
`TwoPhaseStabilityTest.k_values_for_flash`), решает уравнение Рэчфорда-Райза
для мольной доли пара Fv (метод Ньютона с бисекцией как запасным вариантом),
затем итеративно уточняет K-значения методом Ньютона по логарифмам
фугитивностей (аналитический якобиан через `BrusilovskiyEOS._calc_dlogphi_dx_matrix`).
Используется из `VLE.Flash` (двухфазная ветка) и `CompositionalModel`.
"""

import logging
import math

import numpy as np

from calc_core.Composition.Composition import Composition
from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.Utils.EngineConfig import EngineConfig
from calc_core.Utils.Errors import ConvergenceError, InputValidationError
from calc_core.Utils.Validation import validate_positive_pressure, validate_temperature_kelvin

logger = logging.getLogger(__name__)


class PhaseEquilibriumNewton:
    """
    Один расчёт двухфазного равновесия (Rachford-Rice + Ньютон по фугитивностям).

    Экземпляр — одноразовый, под конкретную комбинацию (состав, P, T,
    начальные K). После `find_solve_loop()` держит промежуточные `eos_vapour`/
    `eos_liquid` (объекты `BrusilovskiyEOS` для текущих фазовых составов) как
    атрибуты.
    """

    _RR_EPS = 1e-12
    _RR_NEWTON_TOL = 1e-10
    _RR_MAX_ITER = 1000
    _FUG_MAX_ITER = 1000  # лимит итераций внешнего цикла Ньютона по ln(K)

    def __init__(
        self, composition: Composition, p: float, t: float, k_values,
        *, config: EngineConfig | None = None,
    ):
        """
        Parameters
        ----------
        composition : Composition
            Исходный (фидовый) состав.
        p : float
            Давление, бар.
        t : float
            Температура, K.
        k_values : dict
            Начальное приближение констант равновесия `{компонент: K_i}`,
            обычно из `TwoPhaseStabilityTest.k_values_for_flash`.
        """
        validate_positive_pressure(p, name='p')
        validate_temperature_kelvin(t, name='t')

        if k_values is None:
            raise InputValidationError('Начальные K-значения не должны быть None.')
        try:
            provided_components = set(k_values)
        except TypeError as exc:
            raise InputValidationError(
                'Начальные K-значения должны быть словарём {компонент: K_i}.'
            ) from exc

        expected_components = set(composition.composition)
        if provided_components != expected_components:
            missing = sorted(expected_components - provided_components)
            extra = sorted(provided_components - expected_components)
            raise InputValidationError(
                f'Компоненты K-значений не совпадают с составом: missing={missing}, extra={extra}'
            )

        validated_k = {}
        for component in composition.composition:
            value = k_values[component]
            try:
                number = float(value)
            except (TypeError, ValueError) as exc:
                raise InputValidationError(
                    f'K[{component}] должно быть числом. Передано: {value!r}'
                ) from exc
            if isinstance(value, bool) or not math.isfinite(number) or number <= 0.0:
                raise InputValidationError(
                    f'K[{component}] должно быть конечным и больше 0. Передано: {value!r}'
                )
            validated_k[component] = number

        self._composition = composition
        self._p = float(p)
        self._t = float(t)
        self.config = config or EngineConfig.defaults()
        self._RR_EPS = self.config.flash_rr_epsilon
        self._RR_NEWTON_TOL = self.config.flash_rr_newton_tolerance
        self._RR_MAX_ITER = self.config.flash_rr_max_iterations
        self._FUG_MAX_ITER = self.config.flash_fugacity_max_iterations

        self.zi = composition.composition
        self.k_values = validated_k

        self.L = 0.5
        self.fv = 1.0 - self.L

        self.xi_l = None
        self.yi_v = None
        self.ri = None

        self.eos_vapour = None
        self.eos_liquid = None

        self.convergence = False
        self.trivial_solution = False

        # Внутреннее векторное представление
        self._components = tuple(self.zi.keys())
        self._component_index = {comp: i for i, comp in enumerate(self._components)}
        self._nc = len(self._components)

        self._z_feed = np.fromiter(
            (self.zi[c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )

        self._k = np.fromiter(
            (self.k_values[c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )
        self._log_k = np.log(self._k)

        self._xi_l_arr = None
        self._yi_v_arr = None
        self._ri_arr = None

    # =====================================================================================
    # ВСПОМОГАТЕЛЬНЫЕ МЕТОДЫ
    # =====================================================================================

    def _array_to_dict(self, arr: np.ndarray):
        """Конвертирует numpy-массив (по индексам `self._components`) в `{компонент: значение}`."""
        return {comp: float(arr[i]) for i, comp in enumerate(self._components)}

    def _sync_k_array_from_dict(self):
        """Пересобирает `self._k`/`self._log_k` (numpy) из `self.k_values` (dict) — на случай ручного редактирования dict-версии."""
        self._k = np.fromiter(
            (self.k_values[c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )
        self._log_k = np.log(self._k)

    def _sync_k_dict_from_array(self):
        """Пересобирает `self.k_values` (dict) из `self._k` (numpy) — вызывается после каждого шага Ньютона."""
        self.k_values = {comp: float(self._k[i]) for i, comp in enumerate(self._components)}

    @staticmethod
    def _safe_denominator(arr: np.ndarray, eps: float):
        """
        Защита от деления на (около-)ноль: заменяет элементы `arr`, чьи
        модули меньше `eps`, на `±eps` (знак сохраняется).

        Parameters
        ----------
        arr : np.ndarray
        eps : float
            Порог; значения с `abs(x) < eps` заменяются.

        Returns
        -------
        np.ndarray
            Копия `arr` (исходный массив не модифицируется).
        """
        out = arr.copy()
        mask = np.abs(out) < eps
        if np.any(mask):
            out[mask] = np.where(out[mask] >= 0.0, eps, -eps)
        return out

    # =====================================================================================
    # УРАВНЕНИЕ РЭЧФОРДА-РАЙЗА
    # =====================================================================================

    def _rr_bounds(self):
        """
        Границы допустимых значений Fv для уравнения Рэчфорда-Райза, исходя
        из min/max текущих K-значений: `[1/(1-K_max), 1/(1-K_min)]`.

        Returns
        -------
        tuple[float, float]
            `(fv_min, fv_max)`.

        Raises
        ------
        InputValidationError
            Если K-значения не удовлетворяют требованию `K_min < 1 < K_max`
            (необходимое условие существования двухфазного решения).
        """
        k_min = float(np.min(self._k))
        k_max = float(np.max(self._k))

        if not ((k_min < 1.0) and (k_max > 1.0)):
            raise InputValidationError(
                "Константы равновесия не удовлетворяют требованиям уравнения Рэчфорда-Райза"
            )

        fv_min = 1.0 / (1.0 - k_max)
        fv_max = 1.0 / (1.0 - k_min)
        return fv_min, fv_max

    def _rr_sum(self, fv: float):
        """
        Значение функции Рэчфорда-Райза `sum(z_i (K_i-1) / (1 + Fv(K_i-1)))`
        в точке `fv`. Корень этой функции по Fv — искомое решение.

        Parameters
        ----------
        fv : float

        Returns
        -------
        float
        """
        km1 = self._k - 1.0
        denom = 1.0 + fv * km1
        denom = self._safe_denominator(denom, self._RR_EPS)
        return float(np.sum(self._z_feed * km1 / denom))

    def _rr_sum_and_derivative(self, fv: float):
        """
        Значение функции Рэчфорда-Райза и её производной по Fv в точке `fv`
        (для шага Ньютона в `find_solve_newton`).

        Parameters
        ----------
        fv : float

        Returns
        -------
        tuple[float, float]
            `(значение, производная)`.
        """
        km1 = self._k - 1.0
        denom = 1.0 + fv * km1
        denom = self._safe_denominator(denom, self._RR_EPS)

        rr_sum = np.sum(self._z_feed * km1 / denom)
        rr_der = -np.sum(self._z_feed * (km1 ** 2) / (denom ** 2))
        return float(rr_sum), float(rr_der)

    def find_solve_newton(self):
        """
        Решает уравнение Рэчфорда-Райза методом Ньютона с бисекцией как
        запасным шагом: если шаг Ньютона выходит за текущие границы
        `[fv_left, fv_right]` (или производная слишком мала), делается шаг
        бисекции вместо него. Границы сужаются на каждой итерации по знаку
        `_rr_sum`. Стартовая точка — `self.fv` (по умолчанию 0.5), обрезанная
        в исходные границы `_rr_bounds()`.

        Returns
        -------
        float
            Найденная мольная доля паровой фазы Fv.

        Raises
        ------
        ConvergenceError
            Если `_RR_MAX_ITER` исчерпан без достижения допуска.
        """
        fv_left, fv_right = self._rr_bounds()
        fv = float(np.clip(self.fv, fv_left, fv_right))

        for _ in range(self._RR_MAX_ITER):
            rr_sum, rr_der = self._rr_sum_and_derivative(fv)

            if rr_sum > 0.0:
                fv_left = fv
            else:
                fv_right = fv

            if np.isfinite(rr_der) and abs(rr_der) > self._RR_EPS:
                fv_new = fv - rr_sum / rr_der
            else:
                fv_new = 0.5 * (fv_left + fv_right)

            if (not np.isfinite(fv_new)) or (fv_new <= fv_left) or (fv_new >= fv_right):
                fv_new = 0.5 * (fv_left + fv_right)

            residual_step = abs(fv_new - fv)
            residual_func = abs(self._rr_sum(fv_new))

            fv = fv_new

            if (residual_step <= self._RR_NEWTON_TOL) and (residual_func <= self._RR_NEWTON_TOL):
                return fv

        raise ConvergenceError(
            f'Rachford-Rice (Newton) не сошёлся за {self._RR_MAX_ITER} итераций '
            f'(P={self._p} бар, T={self._t} K, residual={abs(self._rr_sum(fv)):.3e})'
        )

    def find_solve_bisection_v4(self, tol=None):
        """
        Альтернативное (чисто бисекционное, без Ньютона) решение уравнения
        Рэчфорда-Райза. В основном цикле (`find_solve_loop`) не используется —
        там применяется `find_solve_newton`; оставлена как самостоятельный
        запасной метод.

        Parameters
        ----------
        tol : float, optional
            Порог сходимости и по значению функции, и по ширине интервала.
            По умолчанию `Constants.TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE`.

        Returns
        -------
        float
            Найденная мольная доля паровой фазы Fv.

        Raises
        ------
        ConvergenceError
            Если `_RR_MAX_ITER` исчерпан без достижения допуска.
        """
        tol = self.config.flash_bisection_tolerance if tol is None else tol
        fv_left, fv_right = self._rr_bounds()
        f_left = self._rr_sum(fv_left)

        fv = 0.5 * (fv_left + fv_right)

        for _ in range(self._RR_MAX_ITER):
            fv = 0.5 * (fv_left + fv_right)
            f_mid = self._rr_sum(fv)

            if abs(f_mid) < tol or abs(fv_right - fv_left) < tol:
                return fv

            if f_left * f_mid < 0.0:
                fv_right = fv
            else:
                fv_left = fv
                f_left = f_mid

        raise ConvergenceError(
            f'Rachford-Rice (bisection) не сошёлся за {self._RR_MAX_ITER} итераций '
            f'(P={self._p} бар, T={self._t} K, residual={abs(self._rr_sum(fv)):.3e})'
        )

    # =====================================================================================
    # ФАЗОВЫЕ СОСТАВЫ
    # =====================================================================================

    def define_xi_l_yi_v(self):
        """
        По текущим `self.fv`/`self.L` (доля пара/жидкости) и K-значениям
        вычисляет составы фаз: `x_i = z_i / (L + (1-L) K_i)`, `y_i = K_i x_i`.
        Обновляет `self.xi_l`/`self.yi_v` (dict) и их numpy-версии.

        Returns
        -------
        tuple[dict, dict]
            `(xi_l, yi_v)` — составы жидкой и паровой фазы.
        """
        denom = self.L + (1.0 - self.L) * self._k
        denom = self._safe_denominator(denom, self._RR_EPS)

        self._xi_l_arr = self._z_feed / denom
        self._yi_v_arr = self._k * self._xi_l_arr

        self.xi_l = self._array_to_dict(self._xi_l_arr)
        self.yi_v = self._array_to_dict(self._yi_v_arr)

        return self.xi_l, self.yi_v

    # =====================================================================================
    # НЬЮТОН ПО ФУГИТИВНОСТЯМ
    # =====================================================================================

    def fill_jacobian_fug_only(self):
        """
        Аналитический якобиан системы уравнений равновесия фугитивностей
        по ln(K) — правая часть шага Ньютона в `newton_algorithm_fug_only`.

        J_ik = d ln(phi_i^L)/d ln(K_k) - d ln(phi_i^V)/d ln(K_k) - delta_ik

        d ln(phi_i)/d ln(K_k) = sum_j d ln(phi_i)/d x_j * d x_j/d ln(K_k)
        или аналогично для vapor-фазы.

        Здесь:
            dlogphi_dx_l = d ln(phi_i^L) / d x_j^L
            dlogphi_dx_v = d ln(phi_i^V) / d y_j^V
        берутся напрямую из:
            BrusilovskiyEOS._calc_dlogphi_dx_matrix()

        Requires
        --------
        `self.eos_liquid`/`self.eos_vapour` уже посчитаны (`calc_eos()` вызван)
        для текущих фазовых составов — заполняются в `find_solve_loop`.

        Returns
        -------
        np.ndarray, shape (nc, nc)
        """
        # Матрицы производных ln(phi_i) по составам фаз
        dlogphi_dx_l = self.eos_liquid._calc_dlogphi_dx_matrix()   # shape (nc, nc)
        dlogphi_dx_v = self.eos_vapour._calc_dlogphi_dx_matrix()   # shape (nc, nc)

        z_safe = self._safe_denominator(self._z_feed, self._RR_EPS)

        # d x_k / d ln(K_k)
        dx_dlnk = -self._k * (self._xi_l_arr ** 2) * self.fv / z_safe   # shape (nc,)
        # d y_k / d ln(K_k)
        dy_dlnk = self._k * (self._xi_l_arr + dx_dlnk)                  # shape (nc,)

        # Цепное правило:
        # J_liq[i, k] = sum_j dlogphi_dx_l[i, j] * d x_j / d lnK_k
        # Но d x_j / d lnK_k = 0 при j != k, поэтому просто умножение столбцов
        j_liq = dlogphi_dx_l * dx_dlnk[None, :]
        j_vap = dlogphi_dx_v * dy_dlnk[None, :]

        J = j_liq - j_vap - np.eye(self._nc, dtype=np.float64)
        return J

    def fill_column_vector_fug_only(self):
        """
        Вектор невязки `b = ln(phi_i^L) - ln(phi_i^V) - ln(K_i)` — правая
        часть линейной системы `J·delta = -b` шага Ньютона. При сходимости
        (K соответствует равным фугитивностям фаз) стремится к нулю.

        Returns
        -------
        np.ndarray, shape (nc, 1)
        """
        ln_phi_l = self.eos_liquid.get_fugacity_coef_vector_by_root(self.eos_liquid.z)
        ln_phi_v = self.eos_vapour.get_fugacity_coef_vector_by_root(self.eos_vapour.z)

        b = ln_phi_l - ln_phi_v - self._log_k
        return b.reshape(-1, 1)

    def newton_algorithm_fug_only(self):
        """
        Один шаг Ньютона по ln(K): решает `J·delta = -b` (см.
        `fill_jacobian_fug_only`/`fill_column_vector_fug_only`) и обновляет
        `self._log_k`/`self._k`/`self.k_values`.
        """
        J = self.fill_jacobian_fug_only()
        b = self.fill_column_vector_fug_only()

        delta = -np.linalg.solve(J, b).ravel()

        self._log_k = self._log_k + delta
        self._k = np.exp(self._log_k)
        self._sync_k_dict_from_array()

    # =====================================================================================
    # Ri И КРИТЕРИИ ОСТАНОВКИ
    # =====================================================================================

    def calc_Ri(self, eos_vapour, eos_liquid):
        """
        Отношение фугитивностей `r_i = f_i^L / f_i^V = exp(ln f_i^L - ln f_i^V)`
        по компонентам — мера "насколько далеко" текущее K от истинного
        равновесия (в равновесии `r_i = 1` для всех i).

        Parameters
        ----------
        eos_vapour, eos_liquid : BrusilovskiyEOS
            Уже посчитанные (`calc_eos()` вызван) объекты EOS для текущих
            составов фаз.

        Returns
        -------
        dict
            `{компонент: r_i}`.
        """
        ln_f_l = eos_liquid.fugacities   # shape (nc,)
        ln_f_v = eos_vapour.fugacities   # shape (nc,)

        self._ri_arr = np.exp(ln_f_l - ln_f_v)
        self.ri = self._array_to_dict(self._ri_arr)
        return self.ri

    def check_convergence_ri(self, e=None):
        """
        Критерий сходимости по фугитивностям: `sum((r_i - 1)^2) < e`.
        Обновляет и возвращает `self.convergence`.

        Parameters
        ----------
        e : float, optional
            Порог. По умолчанию `Constants.TOL_TWO_PHASE_FLASH_CONVERGENCE`.

        Returns
        -------
        bool
        """
        e = self.config.flash_convergence_tolerance if e is None else e
        sum_ri = float(np.sum((self._ri_arr - 1.0) ** 2))

        if sum_ri < e:
            self.convergence = True
            return True

        self.convergence = False
        return False

    def check_trivial_solution(self):
        """
        Критерий "тривиального решения": `sum(ln(K_i)^2) < tol`, т.е. все
        K_i -> 1 (обе фазы сходятся к одному и тому же составу — на самом
        деле однофазная система). Обновляет и возвращает `self.trivial_solution`.

        Returns
        -------
        bool
        """
        trivial_metric = float(np.sum(self._log_k ** 2))

        if trivial_metric < self.config.flash_trivial_tolerance:
            self.trivial_solution = True
            return True

        self.trivial_solution = False
        return False

    # =====================================================================================
    # ОСНОВНОЙ ЦИКЛ
    # =====================================================================================

    def find_solve_loop(self):
        """
        Главная точка входа: полный цикл решения двухфазного равновесия.

        На каждой итерации: решает Rachford-Rice (`find_solve_newton`) →
        строит составы фаз (`define_xi_l_yi_v`) → считает `BrusilovskiyEOS`
        отдельно для паровой и жидкой фазы → проверяет сходимость по
        фугитивностям (`check_convergence_ri`); если не сошлось — делает шаг
        Ньютона по ln(K) (`newton_algorithm_fug_only`) и проверяет
        тривиальное решение (`check_trivial_solution`). Останов — по любому
        из трёх условий (сходимость / тривиальное решение / лимит 1000
        итераций).

        Returns
        -------
        dict
            `{"yi_v", "xi_l", "Ki", "Fv", "Fl", "Z_v", "Z_l"}` — итоговые
            составы фаз, K-значения, мольные доли пара/жидкости и их
            Z-факторы. Используется напрямую как `phase_equil_result`
            в `VLE.Flash.calculate()`.

        Raises
        ------
        ConvergenceError
            Если внешний цикл Ньютона по ln(K) не сошёлся за `_FUG_MAX_ITER`
            итераций (раньше возвращался последний несошедшийся результат —
            теперь это явный сигнал). Тривиальное решение (K->1) и штатная
            сходимость по фугитивностям — не ошибка, возвращают dict как обычно.
        """
        i = 0

        while True:
            self.fv = self.find_solve_newton()
            self.L = 1.0 - self.fv

            self.define_xi_l_yi_v()

            vapour_composition = self._composition.new_composition(self.yi_v)
            self.eos_vapour = BrusilovskiyEOS(
                composition=vapour_composition,
                p=self._p,
                t=self._t,
            )
            self.eos_vapour.calc_eos()

            liquid_composition = self._composition.new_composition(self.xi_l)
            self.eos_liquid = BrusilovskiyEOS(
                composition=liquid_composition,
                p=self._p,
                t=self._t,
            )
            self.eos_liquid.calc_eos()

            self.calc_Ri(self.eos_vapour, self.eos_liquid)

            self.check_convergence_ri()
            if self.convergence:
                logger.debug("PhaseEquilibriumNewton: сошлось за %d итераций, Fv=%.6f", i, self.fv)
                break

            self.newton_algorithm_fug_only()

            self.check_trivial_solution()
            if self.trivial_solution:
                logger.debug("PhaseEquilibriumNewton: тривиальное решение (K->1) на итерации %d", i)
                break

            i += 1
            if i >= self._FUG_MAX_ITER:
                logger.warning(
                    "PhaseEquilibriumNewton: не сошлось за %d итераций (P=%s, T=%s), последний Fv=%.6f",
                    self._FUG_MAX_ITER, self._p, self._t, self.fv,
                )
                raise ConvergenceError(
                    f"Двухфазное равновесие не сошлось за {self._FUG_MAX_ITER} итераций "
                    f"(P={self._p} бар, T={self._t} K)"
                )

        return {
            "yi_v": self.yi_v,
            "xi_l": self.xi_l,
            "Ki": self.k_values,
            "Fv": self.fv,
            "Fl": self.L,
            "Z_v": self.eos_vapour.z,
            "Z_l": self.eos_liquid.z,
        }
