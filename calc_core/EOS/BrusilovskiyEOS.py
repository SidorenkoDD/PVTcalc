"""
Кубическое уравнение состояния Брусиловского.

Единственная живая реализация EOS в проекте (обобщение Пенга-Робинсона/
Соаве-Редлиха-Квонга, распространено в русской PVT-практике). Используется
везде, где нужны фугитивности/Z-фактор смеси: `PhaseStability.TwoPhaseStabilityTest`
(тест стабильности), `VLE.PhaseEquilibriumNewton` (уточнение K-значений по
фугитивностям), `Utils.FluidPropertiesCalculator` (молярный объём/плотность/Z).

Принимает на вход `Composition` (состав + предвычисленные параметры EOS a/b/c/d
и BIP-матрицу, см. `Composition.evaluate_composition_data`) и условия (P, T).
Отдаёт наружу выбранный Z-фактор (`.z`) и вектор логарифмов фугитивностей
компонент (`.fugacities`) для корня, минимизирующего приведённую энергию
Гиббса — при нескольких действительных корнях кубического уравнения (что
случается вблизи критической точки) это и есть механизм выбора физичного
корня.
"""

import logging
import math
import numpy as np

from calc_core.Composition.Composition import Composition
from calc_core.EOS.BaseEOS import EOS
from calc_core.Utils.Constants import CONSTANT_R
from calc_core.Utils.Cardano import cubic_roots_cardano

logger = logging.getLogger(__name__)


class BrusilovskiyEOS(EOS):
    """
    Векторизованная реализация EOS Брусиловского.
    Внутри класса все вычисления выполняются на numpy-массивах.

    Основные массивы:
        self._x                 : (nc,)   мольные доли смеси
        self._a, _b, _c, _d     : (nc,)   параметры компонент
        self._shift             : (nc,)   shift_parameter
        self._bip               : (nc,nc) матрица бинарных коэффициентов

        self._A, _B, _C, _D     : (nc,)   безразмерные коэффициенты EOS
        self._Aij               : (nc,nc) матрица смешивания

        self.real_roots_eos         : (nr,)
        self.fugacity_by_roots      : (nr,nc)   ln(fi)
        self.fugacity_coef_by_roots : (nr,nc)   ln(phi_i)
        self.normalized_gibbs_energy: (nr,)

    Типичное использование::

        eos = BrusilovskiyEOS(composition, p=100.0, t=350.0)
        z, ln_f = eos.calc_eos()   # z-фактор и ln(fi) выбранного корня
    """

    def __init__(self, composition: Composition, p: float, t: float):
        """
        Parameters
        ----------
        composition : Composition
            Состав с уже посчитанными параметрами EOS (`composition_data['a'|'b'|'c'|'d']`,
            `'bip'`, `'peneloux_correction']`) — то есть после
            `Composition.evaluate_composition_data(...)`.
        p : float
            Давление, бар.
        t : float
            Температура, K.
        """
        super().__init__(composition, p, t)

        self._composition = composition
        self._p = float(p)
        self._t = float(t)

        self.shift_parametr = None
        self.choosen_fugacities = None
        self.choosen_eos_root = None
        self.normalized_gibbs_energy = None
        self.fugacity_by_roots = None
        self.fugacity_coef_by_roots = None
        self.real_roots_eos = None

        # -------------------------------
        # Компоненты и индексы
        # -------------------------------
        self._components = tuple(composition.composition.keys())
        self._component_index = {comp: i for i, comp in enumerate(self._components)}
        self._nc = len(self._components)

        # -------------------------------
        # Состав смеси -> вектор
        # -------------------------------
        self._zi_dict = composition.composition
        self._x = np.fromiter(
            (self._zi_dict[c] for c in self._components),
            dtype=np.float64,
            count=self._nc
        )

        # -------------------------------
        # Свойства компонент -> векторы/матрицы
        # -------------------------------
        comp_data = composition.composition_data

        self._a = np.fromiter(
            (comp_data['a'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc
        )
        self._b = np.fromiter(
            (comp_data['b'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc
        )
        self._c = np.fromiter(
            (comp_data['c'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc
        )
        self._d = np.fromiter(
            (comp_data['d'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc
        )

        self._cpen = np.fromiter(
            (comp_data['peneloux_correction'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc
        )

        self._bip = np.array(
            [
                [comp_data['bip'][ci][cj] for cj in self._components]
                for ci in self._components
            ],
            dtype=np.float64
        )

        # -------------------------------
        # Кэши/рассчитанные коэффициенты
        # -------------------------------
        self._A = None
        self._B = None
        self._C = None
        self._D = None
        self._Aij = None

        self._A_mixture = None
        self._B_mixture = None
        self._C_mixture = None
        self._D_mixture = None

        self._dAm_dx = None
        self._dz_dx = None
        self._dz_dp = None
        self._dlogphi_dx = None
        self._dlogphi_dp = None

        # Выбранный корень Z и его летучести
        self._z_factor = None
        self._fugacities = None

    # =====================================================================================
    # ВСПОМОГАТЕЛЬНЫЕ МЕТОДЫ
    # =====================================================================================

    def _component_to_index(self, component):
        """
        Приводит имя компонента или уже готовый индекс к индексу в `self._components`.

        Parameters
        ----------
        component : str | int
            Имя компонента ('C1', 'nC4' и т.п.) или уже числовой индекс —
            во втором случае возвращается как есть, без проверки диапазона.

        Returns
        -------
        int
        """
        if isinstance(component, int):
            return component
        return self._component_index[component]

    def _root_to_index(self, root: float):
        """
        Находит позицию заданного значения корня в `self.real_roots_eos`.

        Parameters
        ----------
        root : float
            Значение Z-фактора (корня кубического уравнения), сравнение —
            с допуском `rtol=atol=1e-12`.

        Returns
        -------
        int
            Индекс в `self.real_roots_eos`.

        Raises
        ------
        KeyError
            Если такого корня нет среди `self.real_roots_eos`.
        """
        idx = np.where(np.isclose(self.real_roots_eos, root, rtol=1e-12, atol=1e-12))[0]
        if idx.size == 0:
            raise KeyError(f'Корень {root} не найден в real_roots_eos')
        return int(idx[0])

    @property
    def component_names(self):
        """tuple[str, ...]: имена компонент в порядке, соответствующем всем внутренним массивам."""
        return self._components

    @property
    def composition_vector(self):
        """np.ndarray, shape (nc,): вектор мольных долей смеси (тот же порядок, что `component_names`)."""
        return self._x

    # =====================================================================================
    # РАСЧЕТ КОЭФФИЦИЕНТОВ EOS
    # =====================================================================================

    def _calc_component_coefficients(self):
        """
        Вычисляет векторные коэффициенты A_i, B_i, C_i, D_i и матрицу A_ij.
        """
        if self._A is not None:
            return

        rt = CONSTANT_R * self._t
        self._A = self._a * self._p / (rt ** 2)
        self._B = self._b * self._p / rt
        self._C = self._c * self._p / rt
        self._D = self._d * self._p / rt

        sqrtA = np.sqrt(np.outer(self._A, self._A))
        self._Aij = sqrtA * (1.0 - self._bip)

    def _calc_A_mixture(self):
        """
        Смесевой коэффициент A_m = x^T · A_ij · x (квадратичное правило смешения).
        Кэшируется в `self._A_mixture` после первого вызова.

        Returns
        -------
        float
        """
        if self._A_mixture is None:
            self._calc_component_coefficients()
            self._A_mixture = float(self._x @ self._Aij @ self._x)
        return self._A_mixture

    def _calc_B_mixture(self):
        """
        Смесевой коэффициент B_m = x^T · B (линейное правило смешения).
        Кэшируется в `self._B_mixture` после первого вызова.

        Returns
        -------
        float
        """
        if self._B_mixture is None:
            self._calc_component_coefficients()
            self._B_mixture = float(self._x @ self._B)
        return self._B_mixture

    def _calc_C_mixture(self):
        """
        Смесевой коэффициент C_m = x^T · C (линейное правило смешения).
        Кэшируется в `self._C_mixture` после первого вызова.

        Returns
        -------
        float
        """
        if self._C_mixture is None:
            self._calc_component_coefficients()
            self._C_mixture = float(self._x @ self._C)
        return self._C_mixture

    def _calc_D_mixture(self):
        """
        Смесевой коэффициент D_m = x^T · D (линейное правило смешения).
        Кэшируется в `self._D_mixture` после первого вызова.

        Returns
        -------
        float
        """
        if self._D_mixture is None:
            self._calc_component_coefficients()
            self._D_mixture = float(self._x @ self._D)
        return self._D_mixture

    def _calc_roots_eos(self):
        """
        Строит кубическое уравнение Z^3 + E_0 Z^2 + E_1 Z + E_2 = 0 для смеси
        и решает его через `Cardano.cubic_roots_cardano` (только действительные корни).

        Именно здесь могут получиться несколько действительных корней (обычно
        3, иногда практически совпадающих) — это происходит вблизи критической
        точки и разрешается позже в `_choose_eos_root_by_gibbs_energy`.

        Returns
        -------
        np.ndarray
            Действительные корни Z (1-3 штуки в зависимости от режима).
        """
        A_m = self._calc_A_mixture()
        B_m = self._calc_B_mixture()
        C_m = self._calc_C_mixture()
        D_m = self._calc_D_mixture()

        E_0 = C_m + D_m - B_m - 1.0
        E_1 = A_m - B_m * C_m + C_m * D_m - B_m * D_m - D_m - C_m
        E_2 = -(B_m * C_m * D_m + C_m * D_m + A_m * B_m)

        roots = cubic_roots_cardano(1.0, E_0, E_1, E_2, only_real_roots=True)
        return np.array(roots, dtype=np.float64)

    # =====================================================================================
    # РАСЧЕТ ЛЕТУЧЕСТИ И ЭНЕРГИИ ГИББСА
    # =====================================================================================

    def _calc_sum_xi_Aij(self):
        """
        Вектор s_i = sum_j x_j * A_ij
        """
        return self._Aij @ self._x

    def _calc_fugacity_coef_logarithm_vector(self, root: float):
        """
        Возвращает вектор ln(phi_i) для заданного корня Z.

        Parameters
        ----------
        root : float
            Значение Z-фактора (один из `self.real_roots_eos`).

        Returns
        -------
        np.ndarray, shape (nc,)
            Нулевой вектор, если `root <= B_m` (физически недопустимый корень).
        """
        z = float(root)

        A_m = self._A_mixture
        B_m = self._B_mixture
        C_m = self._C_mixture
        D_m = self._D_mixture

        if z <= B_m:
            return np.zeros(self._nc, dtype=np.float64)

        sum_xi_Aij = self._calc_sum_xi_Aij()

        log_term = math.log((z + C_m) / (z + D_m))
        cd_diff = C_m - D_m

        part1 = -math.log(z - B_m)
        part2 = -A_m / cd_diff * (
            (2.0 * sum_xi_Aij) / A_m - (self._C - self._D) / cd_diff
        ) * log_term
        part3 = self._B / (z - B_m)
        part4 = -A_m / cd_diff * (
            self._C / (z + C_m) - self._D / (z + D_m)
        )

        return part1 + part2 + part3 + part4

    def _calc_fugacity_logarithm_vector(self, root: float):
        """
        Возвращает вектор ln(f_i) и ln(phi_i) для заданного корня Z.

        Parameters
        ----------
        root : float
            Значение Z-фактора (один из `self.real_roots_eos`).

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            `(ln_f, ln_phi)`, оба shape (nc,). Оба нулевые, если `root <= B_m`.
        """
        z = float(root)
        B_m = self._B_mixture

        if z <= B_m:
            zeros = np.zeros(self._nc, dtype=np.float64)
            return zeros, zeros

        ln_phi = self._calc_fugacity_coef_logarithm_vector(root)
        ln_f = ln_phi + np.log(self._x * self._p)
        return ln_f, ln_phi

    def _calc_fugacity_by_roots(self):
        """
        Считает ln(f_i) и ln(phi_i) для *каждого* найденного корня Z
        (см. `self.real_roots_eos`), не только для итогового выбранного.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            `(fugacity, fugacity_coef)`, оба shape (nr, nc), где nr —
            число действительных корней.
        """
        nr = len(self.real_roots_eos)
        fugacity = np.zeros((nr, self._nc), dtype=np.float64)
        fugacity_coef = np.zeros((nr, self._nc), dtype=np.float64)

        for ir, root in enumerate(self.real_roots_eos):
            ln_f, ln_phi = self._calc_fugacity_logarithm_vector(root)
            fugacity[ir, :] = ln_f
            fugacity_coef[ir, :] = ln_phi

        return fugacity, fugacity_coef

    def _calc_normalized_gibbs_energy(self):
        """
        Возвращает массив приведённой энергии Гиббса G/(RT) — по одному
        значению на каждый корень из `self.real_roots_eos`. Корень с
        минимальным значением считается термодинамически устойчивым
        (см. `_choose_eos_root_by_gibbs_energy`).

        Физически недопустимые корни (`root <= 0` или `root <= B_m`)
        получают `+inf`, чтобы никогда не быть выбранными минимумом.

        Returns
        -------
        np.ndarray, shape (nr,)
        """
        B_m = self._B_mixture
        gibbs = np.full(len(self.real_roots_eos), np.inf, dtype=np.float64)

        valid = (self.real_roots_eos > 0.0) & (self.real_roots_eos > B_m)
        if np.any(valid):
            gibbs[valid] = np.sum(self.fugacity_by_roots[valid] * self._x[None, :], axis=1)

        return gibbs

    def _choose_eos_root_by_gibbs_energy(self):
        """
        Выбирает физичный корень EOS — тот, что минимизирует приведённую
        энергию Гиббса (см. `_calc_normalized_gibbs_energy`). Это и есть
        механизм разрешения многокорневого случая (несколько действительных
        корней кубического уравнения, типично вблизи критической точки).

        Returns
        -------
        tuple[float, int]
            `(значение выбранного Z, его индекс в self.real_roots_eos)`.
        """
        idx = int(np.argmin(self.normalized_gibbs_energy))
        return float(self.real_roots_eos[idx]), idx

    # =====================================================================================
    # ПРОИЗВОДНЫЕ
    # =====================================================================================

    def _calc_dAm_dx(self):
        """
        dAm/dx_k = 2 * sum_j x_j * A_kj
        Вектор формы (nc,)
        """
        if self._dAm_dx is None:
            self._dAm_dx = 2.0 * (self._Aij @ self._x)
        return self._dAm_dx

    def _calc_dE0_dx(self):
        """
        Вектор dE0/dx_k — производная свободного члена E0 кубического уравнения
        (см. `_calc_roots_eos`) по мольной доле каждого компонента. Используется
        в `_calc_dz_dx` (неявное дифференцирование кубического уравнения по x).

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        return self._C + self._D - self._B

    def _calc_dE1_dx(self):
        """
        Вектор dE1/dx_k — производная коэффициента E1 кубического уравнения
        по мольной доле каждого компонента. Используется в `_calc_dz_dx`.

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        dAm_dx = self._calc_dAm_dx()
        B_m = self._B_mixture
        C_m = self._C_mixture
        D_m = self._D_mixture

        return (
            dAm_dx
            - self._D
            - self._C
            - C_m * self._B
            - B_m * self._C
            + C_m * self._D
            + D_m * self._C
            - B_m * self._D
            - D_m * self._B
        )

    def _calc_dE2_dx(self):
        """
        Вектор dE2/dx_k — производная свободного члена E2 кубического уравнения
        по мольной доле каждого компонента. Используется в `_calc_dz_dx`.

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        dAm_dx = self._calc_dAm_dx()
        A_m = self._A_mixture
        B_m = self._B_mixture
        C_m = self._C_mixture
        D_m = self._D_mixture

        return -(
            self._B * C_m * D_m
            + self._C * B_m * D_m
            + self._D * C_m * B_m
            + self._C * D_m
            + C_m * self._D
            + dAm_dx * B_m
            + A_m * self._B
        )

    def _calc_dz_dx(self):
        """
        Производная выбранного Z-фактора по мольной доле каждого компонента
        (неявное дифференцирование кубического уравнения EOS). Требует, чтобы
        `self._z_factor` уже был выбран (после `calc_eos()`). Кэшируется
        в `self._dz_dx`, сбрасывается в `calc_eos()`/`clear_cache()`.

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        if self._dz_dx is not None:
            return self._dz_dx

        z = self._z_factor
        E_0 = self._C_mixture + self._D_mixture - self._B_mixture - 1.0
        E_1 = (
            self._A_mixture
            - self._B_mixture * self._C_mixture
            + self._C_mixture * self._D_mixture
            - self._B_mixture * self._D_mixture
            - self._D_mixture
            - self._C_mixture
        )

        numerator = (
            self._calc_dE0_dx() * (z ** 2)
            + self._calc_dE1_dx() * z
            + self._calc_dE2_dx()
        )
        denominator = 3.0 * (z ** 2) + 2.0 * E_0 * z + E_1

        self._dz_dx = -(numerator / denominator)
        return self._dz_dx

    def _calc_dz_dp(self):
        """
        Производная выбранного Z-фактора по давлению (неявное дифференцирование
        кубического уравнения EOS). Требует, чтобы `self._z_factor` уже был
        выбран (после `calc_eos()`). Кэшируется в `self._dz_dp`.

        Returns
        -------
        float
        """
        if self._dz_dp is not None:
            return self._dz_dp

        z = self._z_factor
        A_m = self._A_mixture
        B_m = self._B_mixture
        C_m = self._C_mixture
        D_m = self._D_mixture

        E_0 = C_m + D_m - B_m - 1.0
        E_1 = A_m - B_m * C_m + C_m * D_m - B_m * D_m - D_m - C_m

        E_00 = C_m + D_m - B_m
        E_11 = A_m - 2.0 * B_m * C_m + 2.0 * C_m * D_m - 2.0 * B_m * D_m - D_m - C_m
        E_22 = -(3.0 * B_m * C_m * D_m + 2.0 * C_m * D_m + 2.0 * A_m * B_m)

        numerator = E_00 * (z ** 2) + E_11 * z + E_22
        numerator /= self._p
        denominator = 3.0 * (z ** 2) + 2.0 * E_0 * z + E_1

        self._dz_dp = -(numerator / denominator)
        return self._dz_dp

    def _calc_dlogphi_dx_matrix(self):
        """
        Матрица аналитического якобиана d ln(phi_i) / d x_k для выбранного
        корня. Используется в `PhaseEquilibriumNewton.fill_jacobian_fug_only`
        как основной строительный блок якобиана системы уравнений равновесия
        фугитивностей (метод Ньютона по ln(K)).

        Returns
        -------
        np.ndarray, shape (nc, nc)
        """
        if self._dlogphi_dx is not None:
            return self._dlogphi_dx

        z = self._z_factor

        A_m = self._A_mixture
        B_m = self._B_mixture
        C_m = self._C_mixture
        D_m = self._D_mixture

        B_i = self._B[:, None]               # (i,1)
        C_i = self._C[:, None]
        D_i = self._D[:, None]

        B_k = self._B[None, :]               # (1,k)
        C_k = self._C[None, :]
        D_k = self._D[None, :]

        A_ik = self._Aij                     # (i,k)

        dAm_dx = self._calc_dAm_dx()         # (k,)
        sum_xi_Aij = 0.5 * dAm_dx            # (i,) == sum_j x_j A_ij
        sum_xi_Aij_i = sum_xi_Aij[:, None]   # (i,1)

        dz_dx = self._calc_dz_dx()[None, :]  # (1,k)

        cd_diff = C_m - D_m
        log_term = math.log((z + C_m) / (z + D_m))
        part21 = A_m / cd_diff
        part22_i = (2.0 * sum_xi_Aij_i) / A_m - (C_i - D_i) / cd_diff

        deriv_part1 = (dz_dx - B_k) / (z - B_m)

        deriv_part21 = (
            dAm_dx * cd_diff - A_m * (self._C - self._D)
        ) / (cd_diff ** 2)
        deriv_part21 = deriv_part21[None, :]  # (1,k)

        deriv_part22 = (
            2.0 * (A_ik * A_m - dAm_dx[None, :] * sum_xi_Aij_i) / (A_m ** 2)
            + ((C_i - D_i) * (C_k - D_k)) / (cd_diff ** 2)
        )

        deriv_part23 = (
            (dz_dx + C_k) / (z + C_m)
            - (dz_dx + D_k) / (z + D_m)
        )

        deriv_part2 = (
            deriv_part21 * part22_i * log_term
            + part21 * deriv_part22 * log_term
            + part21 * part22_i * deriv_part23
        )

        deriv_part3 = (-B_i * (dz_dx - B_k)) / ((z - B_m) ** 2)

        part42_i = C_i / (z + C_m) - D_i / (z + D_m)

        deriv_part42 = (
            -C_i * (dz_dx + C_k) / ((z + C_m) ** 2)
            + D_i * (dz_dx + D_k) / ((z + D_m) ** 2)
        )

        deriv_part4 = deriv_part21 * part42_i + deriv_part42 * part21

        self._dlogphi_dx = -deriv_part1 - deriv_part2 + deriv_part3 - deriv_part4
        return self._dlogphi_dx

    def _calc_dlogphi_dp_vector(self):
        """
        Вектор d ln(phi_i) / dp для выбранного корня — производная логарифма
        коэффициента фугитивности по давлению.

        ЧИСЛЕННАЯ (центральная разность), не аналитическая — раньше здесь было
        аналитическое выражение (неявное дифференцирование по цепному правилу
        через Z(P)/A_m(P)/B_m(P)/C_m(P)/D_m(P)/s_i(P)), но оно оказалось
        ошибочным: сверка с производной по конечным разностям от уже
        проверенной `_calc_fugacity_coef_logarithm_vector` показала расхождение
        и по знаку, и по величине для большинства компонент (состав KRSNL,
        T=110°C — например C10: численно +1.9e-4, было аналитически −7.5e-3;
        `_calc_dz_dp`, от которой зависела и та формула, при этом сама сошлась
        с численной производной до 9 знаков — ошибка была именно в сборке
        остальных слагаемых). Единственные потребители — `BubblePointCalculator`/
        `DewPointCalculator` (`calc_core/PhaseEnvelope/*Pressure.py`), вне
        основного pipeline `Flash`/`PhaseEquilibriumNewton` (там свой, рабочий
        якобиан — `calc_d_log_phi_i_dxk`, по составу, а не по давлению) —
        поэтому баг не проявлялся раньше.

        Центральная разность строится вокруг ВЫБРАННОГО корня (`self._z_factor`),
        а не вокруг того, что выбрала бы независимая минимизация энергии Гиббса
        при сдвинутом P — вблизи критической точки состава, где действительных
        корней несколько и они близки, это мог бы оказаться другой физический
        корень; здесь нужна производная именно вдоль текущей ветки, поэтому из
        `real_roots_eos` при P±h берётся корень, ближайший к `self._z_factor`.

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        if self._dlogphi_dp is not None:
            return self._dlogphi_dp

        z_anchor = self._z_factor
        h = max(abs(self._p) * 1e-5, 1e-6)

        def ln_phi_at(p_shifted):
            eos = BrusilovskiyEOS(self._composition, p_shifted, self._t)
            eos.calc_eos()
            roots = eos.real_roots_eos
            nearest_root = roots[np.argmin(np.abs(roots - z_anchor))]
            return eos.get_fugacity_coef_vector_by_root(nearest_root)

        ln_phi_plus = ln_phi_at(self._p + h)
        ln_phi_minus = ln_phi_at(self._p - h)

        self._dlogphi_dp = (ln_phi_plus - ln_phi_minus) / (2.0 * h)
        return self._dlogphi_dp

    def calc_d_log_phi_i_dxk(self, component_i, component_k):
        """
        Скалярный доступ к одному элементу якобиана d ln(phi_i)/dx_k
        (см. `_calc_dlogphi_dx_matrix`), по именам компонент.

        Parameters
        ----------
        component_i, component_k : str | int
            Имена компонент (или индексы) — числитель и знаменатель производной.

        Returns
        -------
        float
        """
        i = self._component_to_index(component_i)
        k = self._component_to_index(component_k)
        return float(self._calc_dlogphi_dx_matrix()[i, k])

    def calc_d_log_phi_i_dp(self, component_i):
        """
        Скалярный доступ к одному элементу d ln(phi_i)/dp (см. `_calc_dlogphi_dp_vector`).

        Parameters
        ----------
        component_i : str | int
            Имя компонента (или индекс).

        Returns
        -------
        float
        """
        i = self._component_to_index(component_i)
        return float(self._calc_dlogphi_dp_vector()[i])

    def _calc_dfi_dxk(self, component_i, component_k):
        """
        Производная фугитивности f_i по мольной доле x_k (не логарифма — самой
        фугитивности), через правило произведения: f_i = p * x_i * phi_i.

        Parameters
        ----------
        component_i, component_k : str | int
            Имена компонент (или индексы).

        Returns
        -------
        float
        """
        i = self._component_to_index(component_i)
        k = self._component_to_index(component_k)

        xi = self._x[i]
        log_phi_i = self.fugacity_coef_by_roots[self._chosen_root_index, i]
        d_log_phi_i_dxk = self._calc_dlogphi_dx_matrix()[i, k]

        phi_i = math.exp(log_phi_i)
        return self._p * phi_i * float(i == k) + self._p * xi * phi_i * d_log_phi_i_dxk

    # =====================================================================================
    # ОСНОВНОЙ РАСЧЕТ
    # =====================================================================================

    def calc_eos(self):
        """
        Главная точка входа: полный пересчёт EOS для текущих состава/P/T.

        Порядок: сброс кэша (`clear_cache`) → решение кубического уравнения
        (`_calc_roots_eos`) → фугитивности по каждому корню
        (`_calc_fugacity_by_roots`) → выбор физичного корня по минимуму
        энергии Гиббса (`_choose_eos_root_by_gibbs_energy`). Результат
        сохраняется в `self.z` / `self.fugacities` (через `self._z_factor` /
        `self._fugacities`) и в `self.choosen_eos_root` / `self.choosen_fugacities`.

        Returns
        -------
        tuple[float, np.ndarray]
            (выбранный Z-фактор, вектор ln(f_i))
        """
        self.clear_cache()

        self.real_roots_eos = self._calc_roots_eos()
        self.fugacity_by_roots, self.fugacity_coef_by_roots = self._calc_fugacity_by_roots()
        self.normalized_gibbs_energy = self._calc_normalized_gibbs_energy()

        self.choosen_eos_root, self._chosen_root_index = self._choose_eos_root_by_gibbs_energy()
        self.choosen_fugacities = self.fugacity_by_roots[self._chosen_root_index].copy()
        self.shift_parametr = self._calc_shift_parametr()

        if len(self.real_roots_eos) > 1:
            logger.debug(
                "calc_eos: %d действительных корня %s, выбран Z=%.6f по мин. энергии Гиббса",
                len(self.real_roots_eos), list(self.real_roots_eos), self.choosen_eos_root,
            )

        self._z_factor = self.choosen_eos_root
        self._fugacities = self.choosen_fugacities

        # Производные далее считаются уже для выбранного корня
        self._dz_dx = None
        self._dz_dp = None
        self._dlogphi_dx = None
        self._dlogphi_dp = None

        return self.choosen_eos_root, self.choosen_fugacities

    def clear_cache(self):
        """
        Сбрасывает все закэшированные промежуточные величины (смесевые
        коэффициенты A/B/C/D, производные, найденные корни и фугитивности)
        перед новым пересчётом. Вызывается автоматически в начале `calc_eos()`.
        """
        self._A = None
        self._B = None
        self._C = None
        self._D = None
        self._Aij = None

        self._A_mixture = None
        self._B_mixture = None
        self._C_mixture = None
        self._D_mixture = None

        self._dAm_dx = None
        self._dz_dx = None
        self._dz_dp = None
        self._dlogphi_dx = None
        self._dlogphi_dp = None

        self.real_roots_eos = None
        self.fugacity_by_roots = None
        self.fugacity_coef_by_roots = None
        self.normalized_gibbs_energy = None

        self.choosen_fugacities = None
        self.choosen_eos_root = None
        self._chosen_root_index = None

    def _calc_shift_parametr(self) -> float:
        """
        Смесевой volume-shift параметр (поправка Пенелу) — средневзвешенное
        по мольным долям `self._cpen` (per-компонентный `peneloux_correction`
        из `composition_data`). Используется в `FluidPropertiesCalculator`
        при расчёте молярного объёма/плотности со сдвигом. Имя метода
        (`_parametr`, транслитерация) — устоявшаяся опечатка в проекте, см.
        CLAUDE.md, не переименовывать без отдельного запроса.

        Returns
        -------
        float
        """
        return float(np.sum(self._x * self._cpen))

    # =====================================================================================
    # СВОЙСТВА
    # =====================================================================================

    @property
    def z(self):
        """
        Возвращает выбранный коэффициент сверхсжимаемости.
        """
        return self._z_factor

    @property
    def fugacities(self):
        """
        Возвращает вектор ln(f_i) для выбранного корня.
        """
        return self._fugacities

    # =====================================================================================
    # ДОПОЛНИТЕЛЬНЫЕ ХЕЛПЕРЫ ДЛЯ УДОБСТВА
    # =====================================================================================

    def get_fugacity_vector_by_root(self, root: float):
        """
        Вектор ln(f_i) для конкретного (не обязательно выбранного) корня Z.

        Parameters
        ----------
        root : float
            Значение Z-фактора из `self.real_roots_eos`.

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        idx = self._root_to_index(root)
        return self.fugacity_by_roots[idx]

    def get_fugacity_coef_vector_by_root(self, root: float):
        """
        Вектор ln(phi_i) для конкретного (не обязательно выбранного) корня Z.

        Parameters
        ----------
        root : float
            Значение Z-фактора из `self.real_roots_eos`.

        Returns
        -------
        np.ndarray, shape (nc,)
        """
        idx = self._root_to_index(root)
        return self.fugacity_coef_by_roots[idx]

    def get_fugacity_dict_by_root(self, root: float):
        """
        То же, что `get_fugacity_vector_by_root`, но в виде `{имя_компонента: ln(f_i)}`.
        Необязательный адаптер для старого кода.

        Parameters
        ----------
        root : float

        Returns
        -------
        dict[str, float]
        """
        vec = self.get_fugacity_vector_by_root(root)
        return {c: float(vec[i]) for i, c in enumerate(self._components)}

    def get_fugacity_coef_dict_by_root(self, root: float):
        """
        То же, что `get_fugacity_coef_vector_by_root`, но в виде `{имя_компонента: ln(phi_i)}`.
        Необязательный адаптер для старого кода.

        Parameters
        ----------
        root : float

        Returns
        -------
        dict[str, float]
        """
        vec = self.get_fugacity_coef_vector_by_root(root)
        return {c: float(vec[i]) for i, c in enumerate(self._components)}


if __name__ == '__main__':
    pass