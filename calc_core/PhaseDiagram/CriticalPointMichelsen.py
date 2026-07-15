"""
Критическая точка методом Хайдемана-Халила/Михельсена, в (T,V), а не в (T,P).

Независимая (T,V)-реализация — не трогает и не заменяет ни один из 4
существующих расчётов критической точки в проекте (`CriticalPoint.py`,
`CriticalProperties_den_v.py`, два grid-scan класса в ноутбуке, см. CLAUDE.md).
Ближайший по духу — `CriticalProperties_den_v.py` (тот же тангенциальный
критерий Хайдемана-Халила, та же конструкция матрицы через собственные
значения), но там независимые переменные — T и P; здесь — T и V, как в книге
Michelsen & Mollerup, гл. 9 ("Practical determination of the critical point"):
это не only иначе, это и есть весь смысл метода — по книге, при (T,P)
"quite accurate initial estimates are required", а (T,V) "is much more
well-behaved".

Практическая выгода (T,V)-подхода в этом проекте конкретно: `BrusilovskiyEOS`
основной путь (`calc_eos()`) фиксирует P и решает кубическое уравнение за Z —
именно эта процедура выбора корня (`_choose_eos_root_by_gibbs_energy`)
неустойчива вблизи критической точки (несколько почти совпадающих
действительных корней). Явная (T,V,n)-оценка (`BrusilovskiyEOS.eval_pv_lnphi`/
`eval_at_volume`, добавлена для этого расчёта) вообще не решает кубику — V
задан, P вычисляется напрямую — поэтому не унаследует эту неустойчивость.

Метод (уравнения 46-55 книги, с упрощением Михельсена — одно направление u,
а не полный тензор третьих производных):

1. M_ij = (∂ln(phi_i)/∂n_j) при фиксированных (T,V) — численно, центральной
   разностью (аналитической (T,V)-производной в проекте нет и не было —
   ближайшая аналитическая производная, `calc_d_log_phi_i_dp` при фиксированном
   (T,x), уже документирована как ошибочная, см. коммит про
   PhaseEnvelopeNewton.py, поэтому численная схема — не компромисс, а
   единственный вариант, которому вообще можно доверять здесь).

   Важно: молярная (не полная) формула EOS (`eval_pv_lnphi`) параметризована
   молярным объёмом v, не полным V. Возмущение ОДНОГО мольного числа n_k при
   ФИКСИРОВАННОМ ПОЛНОМ V меняет n_total = Σn_i, а с ним — и молярный объём
   (v = V/n_total), и мольные доли ВСЕХ компонент (x_i = n_i/n_total), не
   только k-той. Наивная центральная разность по x_k при буквально
   зафиксированном v (без пересчёта под новый n_total) считает производную в
   другом ансамбле (T, v=const), а не в нужном (T, V=const) — расхождение
   первого порядка по δ, не мелкая поправка. `_m_matrix` ниже пересчитывает
   v и x заново на каждом пробном n_total.

2. M*_ij = sqrt(z_i*z_j)*M_ij, симметризация, наименьшее собственное значение
   λ1 и его вектор u (`numpy.linalg.eigh` + `argmin` — тот же приём, что
   `CriticalProperties_den_v.py::_compute_b_c`).

3. c = кубический член вдоль u (ур. 55 книги, второе критериальное условие).

   ПЕРВАЯ версия этого модуля переносила `modified_tangent_plane_criterion`
   из `CriticalProperties_den_v.py` "как есть" (тот функционал справедлив в
   (T,P)-ансамбле — там trial-фаза свободно занимает СВОЙ объём при заданном
   P, поправка не нужна). Эмпирически (см. `tests/regression/
   test_critical_point_michelsen.py`, чистый C1) это НЕ работало: c не
   пересекал ноль в физически осмысленной области (T≈Tc), потому что ур. 46
   книги для tpd(T,V,n) в (T,V)-ансамбле требует ДОПОЛНИТЕЛЬНОЕ слагаемое
   `-(V/RT)(P(T,V,n)-P0)` — trial-фаза здесь искусственно занимает ТОТ ЖЕ
   полный объём V, что и референс, при обычно другом давлении, и это
   слагаемое компенсирует разницу. Пропуск этого слагаемого — не мелкая
   поправка, а качественная ошибка (не та функция).

   Текущая версия считает tpd(T,V,n) ПОЛНОСТЬЮ по ур. 46 (см. `_tpd_scalar`)
   как скаляр напрямую (арифметика из P и ln(phi_i), которые уже надёжно
   считает `eval_pv_lnphi`) — без аналитического градиента, поэтому не нужно
   отдельно выводить производную поправочного слагаемого по n (риск которой
   и заставил отказаться от подхода "как в (T,P)-версии"). Третья производная
   вдоль u берётся СТАНДАРТНОЙ 4-точечной центральной разностью
   `[f(2ε)-2f(ε)+2f(-ε)-f(-2ε)]/(2ε³)` по скалярной tpd(s) — не двухшаговым
   трюком (аналитический градиент в одной точке + одна конечная разность),
   которым пользуется (T,P)-версия, но без риска, которым платит трюк.

4. Вложенный цикл: внутренний Newton по T при фиксированном V (λ1 -> 0),
   внешний Newton по V (c -> 0, при каждом шаге использует T, найденную на
   предыдущем внутреннем решении, как стартовое приближение для следующего).

5. Начальные приближения — как в книге: T0 = Σ z_i Tc_i, V0 = 4*b_m(z)
   (мольная коволюм смеси в "сырых" единицах EOS, БЕЗ ×10-конверсии
   `FluidPropertiesCalculator` — см. CLAUDE.md про обнаруженное расхождение
   единиц; здесь используется тот же "сырой" масштаб, что и во всём
   `BrusilovskiyEOS`, тестами `tests/regression/test_eos_volume_explicit.py`
   подтверждено, что `eval_pv_lnphi` на этом масштабе воспроизводит
   существующий (T,P)-путь).
"""

import logging

import numpy as np

from calc_core.Composition.Composition import Composition
from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.Utils.Constants import CONSTANT_R

logger = logging.getLogger(__name__)


class MichelsenCriticalPointCalculator:
    """
    Критическая точка методом Хайдемана-Халила/Михельсена в (T,V).
    См. докстринг модуля за полным описанием алгоритма.
    """

    def __init__(self, composition: Composition):
        self.composition = composition
        self._components = tuple(composition.composition.keys())
        self._nc = len(self._components)

        self.z = np.fromiter(
            (composition.composition[c] for c in self._components),
            dtype=np.float64, count=self._nc,
        )

        comp_data = composition.composition_data
        self._tc = np.fromiter(
            (comp_data['critical_temperature'][c] for c in self._components),
            dtype=np.float64, count=self._nc,
        )

        self.converged = False
        self.iterations_done = 0
        self.result = None

    # =====================================================================================
    # ОЦЕНКА EOS ПРИ ЗАДАННОЙ T (сырые per-компонентные a/b/c/d/bip)
    # =====================================================================================

    def _raw_eos_params_at_t(self, t: float):
        """Выставляет `composition.T = t` (пересчитывает T-зависимые a/BIP —
        см. `Composition.T.setter`) и возвращает свежие сырые массивы
        a/b/c/d/bip в порядке `self._components`."""
        self.composition.T = t
        cd = self.composition.composition_data
        a = np.fromiter((cd['a'][c] for c in self._components), dtype=np.float64, count=self._nc)
        b = np.fromiter((cd['b'][c] for c in self._components), dtype=np.float64, count=self._nc)
        c = np.fromiter((cd['c'][c] for c in self._components), dtype=np.float64, count=self._nc)
        d = np.fromiter((cd['d'][c] for c in self._components), dtype=np.float64, count=self._nc)
        bip = np.array(
            [[cd['bip'][ci][cj] for cj in self._components] for ci in self._components],
            dtype=np.float64,
        )
        return a, b, c, d, bip

    @staticmethod
    def _lnphi_raw_n(a, b, c, d, bip, n: np.ndarray, t: float, v_total: float):
        """
        `ln(phi_i)` для ПРОБНЫХ мольных чисел `n` (не обязаны суммироваться в
        1) при ПОЛНОМ объёме `v_total` (не молярном!) — корректно пересчитывает
        n_total, мольные доли и молярный объём перед вызовом молярной формулы
        `BrusilovskiyEOS.eval_pv_lnphi` (см. докстринг модуля, п.1, за тем,
        почему это нельзя пропустить).

        Returns
        -------
        tuple[float, np.ndarray]
            `(n_total, ln_phi)`.
        """
        n_total = float(np.sum(n))
        x = n / n_total
        v_molar = v_total / n_total
        _, ln_phi = BrusilovskiyEOS.eval_pv_lnphi(a, b, c, d, bip, x, t, v_molar)
        return n_total, ln_phi

    # =====================================================================================
    # M-МАТРИЦА И НАИМЕНЬШЕЕ СОБСТВЕННОЕ ЗНАЧЕНИЕ
    # =====================================================================================

    def _m_matrix(self, a, b, c, d, bip, t: float, v_total: float, delta: float = 1e-6):
        """
        M_ij = ∂ln(phi_i)/∂n_j при фиксированных (T,V) — центральная разность
        по мольному числу n_j (не по мольной доле x_j напрямую — см. докстринг
        модуля). `delta` — малое абсолютное приращение n_j относительно
        n_total=1 (feed уже нормирован), не требует относительного масштабирования.

        Returns
        -------
        np.ndarray, shape (nc, nc)
        """
        nc = self._nc
        m = np.zeros((nc, nc), dtype=np.float64)

        for k in range(nc):
            n_plus = self.z.copy()
            n_plus[k] += delta
            _, lnphi_plus = self._lnphi_raw_n(a, b, c, d, bip, n_plus, t, v_total)

            n_minus = self.z.copy()
            n_minus[k] -= delta
            _, lnphi_minus = self._lnphi_raw_n(a, b, c, d, bip, n_minus, t, v_total)

            m[:, k] = (lnphi_plus - lnphi_minus) / (2.0 * delta)

        return m

    def _lambda1_and_eigvec(self, a, b, c, d, bip, t: float, v_total: float):
        """
        λ1 (наименьшее собственное значение M*) и его собственный вектор u —
        см. докстринг модуля, п.2. Та же конструкция (симметризация + `eigh`
        + `argmin`), что `CriticalProperties_den_v.py::_compute_b_c`.

        Returns
        -------
        tuple[float, np.ndarray]
        """
        m = self._m_matrix(a, b, c, d, bip, t, v_total)
        sqrtz = np.sqrt(self.z)
        m_star = (sqrtz[:, None] * sqrtz[None, :]) * m
        m_sym = 0.5 * (m_star + m_star.T)

        evals, evecs = np.linalg.eigh(m_sym)
        idx = int(np.argmin(evals))
        u = evecs[:, idx].astype(np.float64)
        u /= float(np.linalg.norm(u))
        return float(evals[idx]), u

    def _tpd_scalar(self, a, b, c, d, bip, t: float, v_total: float, y: np.ndarray, lnf_z: np.ndarray, p0: float):
        """
        tpd(T,V,n) по ур. 46 книги, как скаляр, для пробных мольных чисел `y`
        (не обязаны суммироваться в 1) при полном объёме `v_total` — см.
        докстринг модуля, п.3. `lnf_z`/`p0` — уже посчитанные значения в
        референсной точке (n=z), передаются, чтобы не пересчитывать на
        каждый вызов.

        `y` клипается снизу малым положительным числом перед использованием
        (тот же приём `_clip_mole_numbers`, что `CriticalProperties_den_v.py`
        уже применяет в `modified_tangent_plane_criterion`) — без этого для
        компонент с малой мольной долей z_i возмущение `z_i ± eps*u_i*sqrt(z_i)`
        может дать отрицательное пробное число (обычное дело для тяжёлых
        псевдокомпонент с малой долей на многокомпонентном составе) и `log`
        от него — `NaN`.

        `p_y` тоже клипается снизу малым положительным числом — вдали от
        сходимости (T,V) пробного шага вполне может попасть на неустойчивую
        среднюю ветвь кубической EOS (между бинодалью и спинодалью), где P
        от объёма немонотонно и может быть отрицательным — это не баг, а
        известное свойство кубических EOS; клип не даёт `log` упасть в `NaN`
        на промежуточном шаге поиска, не претендуя на физичность именно в
        этой точке (важна лишь непрерывность сигнала для Newton-шага).
        """
        y = np.maximum(y, 1e-30)
        rt = CONSTANT_R * t
        n_total_y, lnphi_y = self._lnphi_raw_n(a, b, c, d, bip, y, t, v_total)
        x_y = y / n_total_y
        p_y, _ = BrusilovskiyEOS.eval_pv_lnphi(a, b, c, d, bip, x_y, t, v_total / n_total_y)
        p_y_safe = max(p_y, 1e-10)
        lnf_y = lnphi_y + np.log(x_y) + np.log(p_y_safe)

        return float(np.sum(y * (lnf_y - lnf_z))) - (v_total / rt) * (p_y - p0)

    def _c_value(self, a, b, c, d, bip, t: float, v_total: float, u: np.ndarray, eps: float = 1e-3):
        """
        Кубический член вдоль `u` (ур. 55 книги, второе критериальное условие)
        — см. докстринг модуля, п.3. Стандартная 4-точечная центральная
        разность третьей производной скалярной tpd(s) вдоль `s`, `n(s) = z +
        s*u*sqrt(z)`.
        """
        z = self.z
        sqrtz = np.sqrt(z)

        p0, lnphi_z = BrusilovskiyEOS.eval_pv_lnphi(a, b, c, d, bip, z, t, v_total)
        lnf_z = lnphi_z + np.log(z) + np.log(max(p0, 1e-10))

        def tpd_at(s: float) -> float:
            y = z + s * u * sqrtz
            return self._tpd_scalar(a, b, c, d, bip, t, v_total, y, lnf_z, p0)

        f_2p = tpd_at(2.0 * eps)
        f_p = tpd_at(eps)
        f_m = tpd_at(-eps)
        f_2m = tpd_at(-2.0 * eps)

        return (f_2p - 2.0 * f_p + 2.0 * f_m - f_2m) / (2.0 * eps ** 3)

    # =====================================================================================
    # ВНУТРЕННИЙ ЦИКЛ (T при фиксированном V)
    # =====================================================================================

    def _solve_inner_t(self, t0: float, v_total: float, tol: float = 1e-9, max_iter: int = 30, fd_t: float = 0.1):
        """
        Newton по T при фиксированном V, драйвит λ1 -> 0 (шаг численный,
        аналитической dλ1/dT нет). Возвращает лучший найденный T, даже если
        не сошёлся за `max_iter` — внешний цикл сам проверяет `abs(lam1)`.

        Returns
        -------
        tuple[float, np.ndarray, tuple, float]
            `(T, u, (a,b,c,d,bip) при этой T, λ1)`.
        """
        t = float(t0)
        params = None
        lam1 = None
        u = None

        for _ in range(max_iter):
            params = self._raw_eos_params_at_t(t)
            lam1, u = self._lambda1_and_eigvec(*params, t, v_total)

            if abs(lam1) < tol:
                return t, u, params, lam1

            params_p = self._raw_eos_params_at_t(t + fd_t)
            lam1_p, _ = self._lambda1_and_eigvec(*params_p, t + fd_t, v_total)
            params_m = self._raw_eos_params_at_t(t - fd_t)
            lam1_m, _ = self._lambda1_and_eigvec(*params_m, t - fd_t, v_total)

            dlam_dt = (lam1_p - lam1_m) / (2.0 * fd_t)
            if abs(dlam_dt) < 1e-14:
                logger.warning("Внутренний цикл (T): dλ1/dT слишком мал при T=%.4f K, V=%.6f — останов", t, v_total)
                break

            step = -lam1 / dlam_dt
            max_step = 0.25 * t
            if abs(step) > max_step:
                step = max_step * np.sign(step)

            t_new = t + step
            if t_new <= 0:
                t_new = t * 0.5
            t = t_new

        params = self._raw_eos_params_at_t(t)
        lam1, u = self._lambda1_and_eigvec(*params, t, v_total)
        return t, u, params, lam1

    # =====================================================================================
    # ПУБЛИЧНЫЙ МЕТОД (внешний цикл по V)
    # =====================================================================================

    def calculate(
        self,
        t0: float | None = None,
        v0: float | None = None,
        tol_lambda1: float = 1e-9,
        tol_c: float = 1e-7,
        max_outer: int = 40,
        max_inner: int = 30,
        damping: float = 1.0,
    ) -> dict:
        """
        Returns
        -------
        dict
            `{'T': K, 'P': бар, 'V': молярный объём (сырые ед. EOS), 'T_C': °C}`.

        Raises
        ------
        RuntimeError
            Если не сошёлся за `max_outer` итераций даже с демпфированием
            (тот же паттерн отката, что `CriticalProperties_den_v.py::calculate`
            — уменьшение `damping` вдвое и повторная попытка).
        """
        if t0 is None:
            t0 = float(np.dot(self.z, self._tc))
        _, b_raw, _, _, _ = self._raw_eos_params_at_t(t0)
        # `b` в этой EOS-семье не зависит от T (только `a` — см. `BRS_EOS_DB.
        # calc_a_b_c_d`), поэтому b_m(z) — константа на весь расчёт, а не
        # только начальное приближение: P(T,V,z) имеет полюс при V=b_m(z)
        # (см. `eval_pv_lnphi`), и внешний цикл не должен пускать V туда ни
        # на одном шаге — не только в качестве стартовой точки.
        b_m_floor = float(np.dot(self.z, b_raw))
        if v0 is None:
            v0 = 4.0 * b_m_floor

        attempt_damping = damping
        while attempt_damping > 0.01:
            try:
                return self._calculate_once(
                    t0, v0, tol_lambda1, tol_c, max_outer, max_inner, attempt_damping, b_m_floor,
                )
            except RuntimeError as exc:
                logger.warning("Michelsen: попытка с damping=%.3f не сошлась (%s), снижаем", attempt_damping, exc)
                attempt_damping *= 0.5

        raise RuntimeError(
            "MichelsenCriticalPointCalculator не сошёлся ни с одним damping >= 0.01"
        )

    def _calculate_once(self, t0, v0, tol_lambda1, tol_c, max_outer, max_inner, damping, b_m_floor) -> dict:
        v = float(v0)
        t_guess = float(t0)
        v_min = 1.05 * b_m_floor

        for outer_iter in range(max_outer):
            t_conv, u, params, lam1 = self._solve_inner_t(t_guess, v, tol=tol_lambda1, max_iter=max_inner)
            c_val = self._c_value(*params, t_conv, v, u)

            logger.debug(
                "Michelsen outer #%d: T=%.4f K, V=%.6f, λ1=%.3e, c=%.3e",
                outer_iter, t_conv, v, lam1, c_val,
            )

            if abs(lam1) < tol_lambda1 and abs(c_val) < tol_c:
                a, b, c, d, bip = params
                p_final, _ = BrusilovskiyEOS.eval_pv_lnphi(a, b, c, d, bip, self.z, t_conv, v)
                self.converged = True
                self.iterations_done = outer_iter + 1
                self.result = {'T': t_conv, 'P': p_final, 'V': v, 'T_C': t_conv - 273.15}
                logger.info(
                    "Критическая точка (Михельсен): T=%.4f K (%.2f °C), P=%.4f бар, V=%.6f",
                    t_conv, t_conv - 273.15, p_final, v,
                )
                return self.result

            dv = max(v * 1e-4, 1e-8)
            c_plus = self._c_value(*params, t_conv, v + dv, u)
            c_minus = self._c_value(*params, t_conv, v - dv, u)
            dc_dv = (c_plus - c_minus) / (2.0 * dv)

            if abs(dc_dv) < 1e-20:
                raise RuntimeError(f"dc/dV слишком мал при V={v:.6f} (outer_iter={outer_iter})")

            step = -c_val / dc_dv
            max_step = 0.25 * v
            if abs(step) > max_step:
                step = max_step * np.sign(step)

            v_new = v + damping * step
            if v_new <= v_min:
                # Полпути к текущей границе (не жёстко в v_min) — не даёт
                # шагу "залипнуть" ровно на границе на нескольких подряд
                # итерациях (см. докстринг v_min выше — полюс P(T,V,z) при
                # V=b_m(z), сюда идти нельзя ни на шаг ближе).
                v_new = (v + v_min) / 2.0

            v = v_new
            t_guess = t_conv

        raise RuntimeError(f"Достигнуто максимальное число внешних итераций: {max_outer}")
