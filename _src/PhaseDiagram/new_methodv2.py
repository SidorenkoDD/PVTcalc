import numpy as np
import logging
import pandas as pd
import matplotlib.pyplot as plt
from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')

class SaturationPressure:
    def __init__(self, composition: Composition, t_K: float, p_max_bar: float = 1000, p_min_bar: float = 0.01,
                 p_guess_sat: float = None, p_guess_dew: float = None):
        self.composition = composition
        self.composition.T = t_K
        self.t_K = t_K
        self.p_max_bar_sat = p_max_bar
        self.p_min_bar_sat = p_min_bar
        self.p_i = (p_max_bar + p_min_bar) / 2.0

        self.p_min_bar_dew = 0.0
        self.p_max_bar_dew = None
        self.p_i_dew = None
        self._p_guess_dew = p_guess_dew  # Сохраняем предсказание для dew point

        self._last_sum_y_dp  = 0.0
        self._last_Ykz_dp = 0.0
        self.Sum_dp = 0.0

        self.zi = composition.composition
        self._components = tuple(self.zi.keys())
        self._z_arr = np.array(list(self.zi.values()))
        self._nc = len(self._components)
        
        self._last_sum_y = 0.0
        self._last_Ykz = 0.0
        
        # 🔑 ПРИМЕНЯЕМ НАЧАЛЬНЫЕ ПРИБЛИЖЕНИЯ
        if p_guess_sat is not None and p_min_bar < p_guess_sat < p_max_bar:
            self.p_i = p_guess_sat
            self.p_min_bar_sat = max(p_min_bar, p_guess_sat * 0.6)
            self.p_max_bar_sat = min(p_max_bar, p_guess_sat * 1.8)
            logger.info(f"💡 Начальное приближение Psat={p_guess_sat:.3f} bar")
            
        if p_guess_dew is not None and 0 < p_guess_dew < p_max_bar:
            self.p_i_dew = p_guess_dew
            self.p_min_bar_dew = max(0.01, p_guess_dew * 0.6)
            logger.info(f"💡 Начальное приближение Pdew={p_guess_dew:.3f} bar")

        logger.info(f"Инициализация: T={self.t_K:.2f} K, P_min={self.p_min_bar_sat:.4f}, P_max={self.p_max_bar_sat:.4f}")

    def _find_S(self, p: float) -> dict:
        phase_stability = TwoPhaseStabilityTest(self.composition, p, t=self.t_K)
        phase_stability.calculate_phase_stability()
        
        S_l, S_v = phase_stability.S_l, phase_stability.S_v
        
        if (S_l - 1.0)  < 1e-5 and (S_v - 1.0)  < 1e-5:
            return {'S': 0.0, 'yi': np.zeros(self._nc), 'letuch_sp': None, 'letuch_z': None, 'stable': phase_stability.stable}

        yi = np.zeros(self._nc)
        letuch_z = phase_stability._mixture_fugacities_arr.copy()
        letuch_sp = None

        if S_l  > 1.0:
            if S_l  > S_v:
                letuch_sp = phase_stability.liquid_eos.fugacities.copy()
                k_values = phase_stability._k_l_arr.copy()
                yi = self._z_arr / k_values
            else:
                letuch_sp = phase_stability.vapour_eos.fugacities.copy()
                k_values = phase_stability._k_v_arr.copy()
                yi = self._z_arr * k_values
        else:
            if S_v  < 1.0:
                return {'S': 0.0, 'yi': np.zeros(self._nc), 'letuch_sp': None, 'letuch_z': None, 'stable': phase_stability.stable}

        if S_v  > 1.0:
            if S_v  > S_l:
                letuch_sp = phase_stability.vapour_eos.fugacities.copy()
                k_values = phase_stability._k_v_arr.copy()
                yi = self._z_arr * k_values
            elif S_l  < 1.0:
                return {'S': 0.0, 'yi': np.zeros(self._nc), 'letuch_sp': None, 'letuch_z': None, 'stable': phase_stability.stable}

        S_sp = float(np.sum(yi))
        return {'S':  S_sp, 'yi': yi, 'letuch_sp': letuch_sp, 'letuch_z': letuch_z, 'stable': phase_stability.stable}

    def _find_stability_bracket(self, p_low_init=0.01, p_high_init=1000.0) -> tuple[float, float]:
        logger.info(f"🔍 Поиск bracket: start=[{p_low_init:.3f}, {p_high_init:.3f}]")
        res_low = self._find_S(p_low_init)
        res_high = self._find_S(p_high_init)
        stable_low = res_low['stable']
        stable_high = res_high['stable']
        
        if stable_low == stable_high:
            if stable_high:
                p_search = p_high_init
                for step in range(10):
                    p_search *= 0.5
                    if p_search < p_low_init * 0.1:
                        logger.warning("Не найдена unstable точка при снижении")
                        return p_low_init, p_high_init
                    if not self._find_S(p_search)['stable']:
                        p_unstable = p_search
                        p_stable = p_search * 2.0
                        break
                else:
                    logger.warning("Превышен лимит поиска unstable")
                    return p_low_init, p_high_init
            else:
                p_search = p_high_init
                for step in range(10):
                    p_search *= 1.5
                    if p_search > 3000:
                        logger.warning("Не найдена stable точка при повышении")
                        return p_low_init, p_high_init
                    if self._find_S(p_search)['stable']:
                        p_stable = p_search
                        p_unstable = p_search / 1.5
                        break
                else:
                    logger.warning("Превышен лимит поиска stable")
                    return p_low_init, p_high_init
        else:
            if stable_low:
                p_stable, p_unstable = p_low_init, p_high_init
            else:
                p_stable, p_unstable = p_high_init, p_low_init

        if p_stable < p_unstable:
            p_stable, p_unstable = p_unstable, p_stable

        if abs(p_stable - p_unstable) < 1e-6:
            p_unstable *= 0.99
            p_stable *= 1.01

        for i in range(10):
            p_mid = (p_stable + p_unstable) / 2.0
            if self._find_S(p_mid)['stable']:
                p_stable = p_mid
            else:
                p_unstable = p_mid
            if abs(p_stable - p_unstable) < 1e-8:
                break
                
        logger.info(f"Диапазон поиска: unstable={p_unstable:.4f} bar -> stable={p_stable:.4f} bar")
        return p_unstable, p_stable

    def _update_dew_bounds(self):
        """Вспомогательный метод для установки границ Dew Point с учётом предсказания"""
        self.p_max_bar_dew = self.p_i
        if self._p_guess_dew is None:
            self.p_i_dew = self.p_i / 2
        else:
            self.p_i_dew = self._p_guess_dew
            self.p_min_bar_dew = min(self.p_min_bar_dew, self.p_i_dew * 0.5)

    def sp_process(self, lambd: float = 1.0) -> bool | None:
        cur_s_sp = self._find_S(self.p_i)
        while np.isclose(cur_s_sp['S'], 0.0, atol=1e-6):
            self.p_max_bar_sat = self.p_i
            self.p_i = (self.p_max_bar_sat + self.p_min_bar_sat) / 2.0
            if (self.p_max_bar_sat - self.p_min_bar_sat) < 1e-6:
                return None
            cur_s_sp = self._find_S(self.p_i)

        s_sp_val = cur_s_sp['S']
        f_z = cur_s_sp['letuch_z']
        f_sp = cur_s_sp['letuch_sp']
        
        r_sp = np.exp(f_z - f_sp) / s_sp_val
        y_sp = cur_s_sp['yi'] * np.power(r_sp, lambd)
        self._last_sum_y = float(np.sum(y_sp))
        
        ratio = np.clip(y_sp / self._z_arr, 1e-12, None)
        r_sp_safe = np.clip(r_sp, 1e-12, None)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio = np.log(ratio)
            safe_denom = np.where(np.abs(log_ratio) > 1e-12, log_ratio, 1.0)
            self.Sum = float(np.sum(np.log(r_sp_safe) / safe_denom))

        self._last_Ykz = float(np.sum(ratio))
        
        if (abs(1.0 - self._last_sum_y) < 1e-5) or (self._last_Ykz ** 2 < 1e-5):
            return True
        else:
            if self._last_sum_y > 1.0:
                self.p_min_bar_sat = self.p_i
            else:
                self.p_max_bar_sat = self.p_i
            self.p_i = (self.p_max_bar_sat + self.p_min_bar_sat) / 2.0
            return False

    def sp_convergence_loop(self, max_iter: int = 100) -> float | None:
        logger.info("Запуск поиска давления насыщения")
        
        bracket = self._find_stability_bracket(self.p_min_bar_sat, self.p_max_bar_sat)
        if bracket and bracket[0] is not None and abs(bracket[1] - bracket[0]) > 1e-6:
            p_unstable, p_stable = bracket
            self.p_min_bar_sat = p_unstable
            self.p_max_bar_sat = p_stable
            self.p_i = (p_stable + p_unstable) / 2.0
        else:
            logger.warning("Границы не найдены или нулевые, используем исходные границы")
            self.p_min_bar_sat = 0.01
            self.p_max_bar_sat = 1000.0
            self.p_i = 500.0
            
        prev_p = None
        for i in range(1, max_iter + 1):
            if prev_p is not None and abs(self.p_i - prev_p) < 1e-6:
                logger.warning(f"Давление застряло на итерации {i}. P={self.p_i:.6f}")
                break
            prev_p = self.p_i
            
            if (self.p_max_bar_sat - self.p_min_bar_sat) < 1e-6:
                logger.info(f"Интервал схлопнулся на итерации {i}. Возврат P={self.p_i:.4f} bar")
                self._update_dew_bounds()
                return self.p_i
                
            res = self.sp_process()
            if res is None:
                logger.warning(f"sp_process вернул None на итерации {i}")
                self._update_dew_bounds()
                return self.p_i
            if res:
                logger.info(f"Сходимость достигнута на итерации {i}. P_sat={self.p_i:.4f} bar")
                self._update_dew_bounds()
                return self.p_i
                
        logger.warning(f"Превышен лимит итераций ({max_iter}). Последнее P={self.p_i:.4f} bar")
        self._update_dew_bounds()
        return self.p_i

    def dp_process(self, lambd: float = 1.0) -> bool | None:
        cur_s_dp = self._find_S(self.p_i_dew)
        while cur_s_dp['S'] == 0.0:
            self.p_min_bar_dew = self.p_i_dew
            self.p_i_dew = (self.p_max_bar_dew + self.p_min_bar_dew) / 2.0
            if (self.p_max_bar_dew - self.p_min_bar_dew) < 1e-3:
                return None
            cur_s_dp = self._find_S(self.p_i_dew)

        s_dp_val = cur_s_dp['S']
        f_z = cur_s_dp['letuch_z']
        f_dp = cur_s_dp['letuch_sp']
        
        r_dp = np.exp(f_z - f_dp) / s_dp_val
        y_dp = cur_s_dp['yi'] * np.power(r_dp, lambd)
        self._last_sum_y_dp = float(np.sum(y_dp))
        
        ratio = np.clip(y_dp / self._z_arr, 1e-12, None)
        r_dp_safe = np.clip(y_dp, 1e-12, None)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio = np.log(ratio)
            safe_denom = np.where(np.abs(log_ratio) > 1e-12, log_ratio, 1.0)
            self.Sum = float(np.sum(np.log(r_dp_safe) / safe_denom))

        self._last_Ykz_dp = float(np.sum(ratio))

        if (abs(1.0 - self._last_sum_y_dp) < 1e-4) or (self._last_Ykz_dp ** 2 < 1e-4):
            return True
        else:
            self.p_max_bar_dew = self.p_i_dew
            self.p_i_dew = (self.p_max_bar_dew + self.p_min_bar_dew) / 2.0
            return False

    def dp_convergence_loop(self, max_iter: int = 100) -> float | None:
        logger.info("Запуск поиска давления начала конденсации (Dew Point)")
        self.dp_process(lambd=1.0)
        if (self.p_max_bar_dew - self.p_min_bar_dew) < 1e-3:
            return None

        for _ in range(max_iter):
            if (abs(1.0 - self._last_sum_y_dp) < 1e-4) or (self._last_Ykz_dp ** 2 < 1e-4):
                break
            self.dp_process(lambd=1.0)
            if (self.p_max_bar_dew - self.p_min_bar_dew) < 1e-3:
                return None

        self.p_dew = self.p_i_dew
        logger.info(f"Dew Point найден: {self.p_dew:.4f} bar")
        return self.p_dew


class PhaseDiagram:
    def __init__(self, composition: Composition, p_max_bar: float, t_min_C: float, t_max_C: float, t_step_C: float):
        self.composition = composition
        self.p_max_bar = p_max_bar
        self.temps_C = np.arange(t_min_C, t_max_C + t_step_C/2, t_step_C)
        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_bar': []}

    def calculate(self):
        logger.info("Начало расчёта фазовой диаграммы...")
        prev_pb = None
        prev_pdew = None
        
        for t_C in self.temps_C:
            t_K = t_C + 273.15
            
            # 🔑 Передаём предыдущие значения как начальные приближения
            sat_calc = SaturationPressure(
                self.composition, t_K, 
                p_max_bar=self.p_max_bar,
                p_guess_sat=prev_pb,
                p_guess_dew=prev_pdew
            )
            
            pb = sat_calc.sp_convergence_loop()
            pdew = sat_calc.dp_convergence_loop()
            
            self.results['Temp_C'].append(t_C)
            self.results['Bubble_bar'].append(pb if pb is not None else np.nan)
            self.results['Dew_bar'].append(pdew if pdew is not None else np.nan)
            
            # Обновляем предсказания для следующей итерации
            if pb is not None: prev_pb = pb
            if pdew is not None: prev_pdew = pdew

    def plot(self):
        df = pd.DataFrame(self.results)
        if df.empty:
            logger.warning("Нет данных для построения.")
            return

        plt.figure(figsize=(10, 6))
        plt.plot(df['Temp_C'], df['Bubble_bar'], 'b-o', label='Bubble Point', markersize=4, linewidth=2)
        plt.plot(df['Temp_C'], df['Dew_bar'], 'r-s', label='Dew Point', markersize=4, linewidth=2)
        
        valid = df.dropna(subset=['Bubble_bar', 'Dew_bar'])
        if not valid.empty:
            crit_idx = np.abs(valid['Bubble_bar'] - valid['Dew_bar']).idxmin()
            crit_t, crit_p = valid.loc[crit_idx, 'Temp_C'], valid.loc[crit_idx, 'Bubble_bar']
            plt.scatter([crit_t], [crit_p], c='gold', s=150, zorder=5, edgecolor='black', 
                        label=f'Critical ({crit_t:.1f}°C, {crit_p:.1f} bar)')

        plt.xlabel('Temperature (°C)', fontsize=12)
        plt.ylabel('Pressure (bar)', fontsize=12)
        plt.title('Phase Envelope (P-T Diagram)', fontsize=14)
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.show()

    def get_data(self) -> pd.DataFrame:
        return pd.DataFrame(self.results)