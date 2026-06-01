import numpy as np
import logging
import matplotlib.pyplot as plt
import pandas as pd
from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest

logger = logging.getLogger(__name__)

# ==============================================================================
# 1. TOCHKA KIPENIYA (BUBBLE POINT)
# ==============================================================================
class SaturationPressure:
    def __init__(self, composition: Composition, t_K: float, p_max_bar: float = 1000, p_min_bar: float = 0.01):
        self.composition_object = composition
        self.t_K = t_K
        self.p_max_bar_sat = p_max_bar
        self.p_min_bar_sat = p_min_bar
        self.p_i = (p_max_bar + p_min_bar) / 2.0
        
        self.zi = composition.composition
        self._components = tuple(self.zi.keys())
        self._z_arr = np.array(list(self.zi.values()))
        self._nc = len(self._components)
        
        self._last_sum_y = 0.0
        self._last_Ykz = 0.0
        logger.info(f"[Bubble] Инициализация: T={self.t_K:.2f} K, P_min={self.p_min_bar_sat:.4f}, P_max={self.p_max_bar_sat:.4f}")

    def _find_S(self, p: float) -> dict:
        phase_stability = TwoPhaseStabilityTest(self.composition_object, p, t=self.t_K)
        phase_stability.calculate_phase_stability()
        S_l, S_v = phase_stability.S_l, phase_stability.S_v
        
        if (S_l - 1.0) < 1e-5 and (S_v - 1.0) < 1e-5:
            return {'S': 0.0, 'yi': np.zeros(self._nc), 'letuch_sp': None, 'letuch_z': None, 'stable': phase_stability.stable}

        yi = np.zeros(self._nc)
        letuch_z = phase_stability._mixture_fugacities_arr.copy()
        letuch_sp = None

        if S_l > 1.0:
            if S_l > S_v:
                letuch_sp = phase_stability.liquid_eos.fugacities.copy()
                k_values = phase_stability._k_l_arr.copy()
                yi = self._z_arr / k_values
            else:
                letuch_sp = phase_stability.vapour_eos.fugacities.copy()
                k_values = phase_stability._k_v_arr.copy()
                yi = self._z_arr * k_values
        else:
            if S_v < 1.0:
                return {'S': 0.0, 'yi': np.zeros(self._nc), 'letuch_sp': None, 'letuch_z': None, 'stable': phase_stability.stable}

        if S_v > 1.0:
            if S_v > S_l:
                letuch_sp = phase_stability.vapour_eos.fugacities.copy()
                k_values = phase_stability._k_v_arr.copy()
                yi = self._z_arr * k_values
            elif S_l < 1.0:
                return {'S': 0.0, 'yi': np.zeros(self._nc), 'letuch_sp': None, 'letuch_z': None, 'stable': phase_stability.stable}

        return {'S': float(np.sum(yi)), 'yi': yi, 'letuch_sp': letuch_sp, 'letuch_z': letuch_z, 'stable': phase_stability.stable}

    def _find_stability_bracket(self, p_low_init=0.01, p_high_init=1000.0) -> tuple[float, float]:
        logger.info(f"[Bubble] 🔍 Сужение bracket: start=[{p_low_init:.3f}, {p_high_init:.3f}]")
        p_high = p_high_init
        last_stable_p = p_high
        
        for _ in range(20):
            res = self._find_S(p_high)
            if not res['stable']:
                p_unstable = p_high
                p_stable = last_stable_p
                break
            last_stable_p = p_high
            p_high *= 0.5
        else:
            p_high = p_high_init
            for _ in range(15):
                p_high *= 1.5
                if p_high > 2500: break
                if self._find_S(p_high)['stable']:
                    p_stable, p_unstable = p_high, p_high_init
                    break
            else:
                logger.warning("[Bubble] ⚠️ Переход stable/unstable не найден. Использую исходные границы.")
                return p_low_init, p_high_init

        for _ in range(8):
            p_mid = (p_stable + p_unstable) / 2.0
            if self._find_S(p_mid)['stable']:
                p_stable = p_mid
            else:
                p_unstable = p_mid
                
        logger.info(f"[Bubble] ✅ Bracket: unstable={p_unstable:.4f} -> stable={p_stable:.4f}")
        return p_unstable, p_stable

    def sp_process(self, lambd: float = 1.0) -> bool | None:
        cur_s_sp = self._find_S(self.p_i)
        while np.isclose(cur_s_sp['S'], 0.0, atol=1e-6):
            self.p_max_bar_sat = self.p_i
            self.p_i = (self.p_max_bar_sat + self.p_min_bar_sat) / 2.0
            if (self.p_max_bar_sat - self.p_min_bar_sat) < 1e-12: return None
            cur_s_sp = self._find_S(self.p_i)

        s_sp_val = cur_s_sp['S']
        f_z, f_sp = cur_s_sp['letuch_z'], cur_s_sp['letuch_sp']
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
        logger.debug(f"[Bubble] P={self.p_i:.4f} | sum_y={self._last_sum_y:.6f} | Ykz²={self._last_Ykz**2:.2e}")

        if (abs(1.0 - self._last_sum_y) < 1e-4) or (self._last_Ykz ** 2 < 1e-4):
            return True
        else:
            self.p_min_bar_sat = self.p_i
            self.p_i = (self.p_max_bar_sat + self.p_min_bar_sat) / 2.0
            return False

    def sp_convergence_loop(self, max_iter: int = 100) -> float | None:
        bracket = self._find_stability_bracket(self.p_min_bar_sat, self.p_max_bar_sat)
        if bracket[0] is not None:
            p_unstable, p_stable = bracket
            self.p_min_bar_sat = p_unstable
            self.p_max_bar_sat = p_stable
            self.p_i = (p_stable + p_unstable) / 2.0

        for i in range(1, max_iter + 1):
            if (self.p_max_bar_sat - self.p_min_bar_sat) < 1e-4:
                return self.p_i
            res = self.sp_process()
            if res is None or res:
                return self.p_i
        logger.warning(f"[Bubble] ⚠️ Лимит итераций. P={self.p_i:.4f}")
        return self.p_i


# ==============================================================================
# 2. TOCHKA ROSY (DEW POINT)
# ==============================================================================
class DewPointPressure:
    def __init__(self, composition: Composition, t_K: float, p_max_bar: float = 1000, p_min_bar: float = 0.01):
        self.composition_object = composition
        self.t_K = t_K
        self.p_max_dew = p_max_bar
        self.p_min_dew = p_min_bar
        self.p_i = (p_max_bar + p_min_bar) / 2.0
        
        self.zi = composition.composition
        self._components = tuple(self.zi.keys())
        self._z_arr = np.array(list(self.zi.values()))
        self._nc = len(self._components)
        
        self._last_sum_y = 0.0
        self._last_Ykz = 0.0
        logger.info(f"[Dew] Инициализация: T={self.t_K:.2f} K, P_min={self.p_min_dew:.4f}, P_max={self.p_max_dew:.4f}")

    def _find_S(self, p: float) -> dict:
        # Физика теста стабильности идентична для bubble/dew
        return SaturationPressure._find_S(SaturationPressure(self.composition_object, self.t_K, self.p_max_dew, self.p_min_dew), p)

    def _find_stability_bracket(self, p_low_init=0.01, p_high_init=1000.0) -> tuple[float, float]:
        logger.info(f"[Dew] 🔍 Сужение bracket: start=[{p_low_init:.3f}, {p_high_init:.3f}]")
        p_high = p_high_init
        last_stable_p = p_high
        
        for step in range(20):
            res = self._find_S(p_high)
            if not res['stable']:
                p_unstable = p_high
                p_stable = last_stable_p
                break
            last_stable_p = p_high
            p_high *= 0.5
        else:
            p_high = p_high_init
            for step in range(15):
                p_high *= 1.5
                if p_high > 2500: break
                if self._find_S(p_high)['stable']:
                    p_stable, p_unstable = p_high, p_high_init
                    break
            else:
                logger.warning("[Dew] ⚠️ Переход stable/unstable не найден.")
                return p_low_init, p_high_init

        for _ in range(8):
            p_mid = (p_stable + p_unstable) / 2.0
            if self._find_S(p_mid)['stable']:
                p_stable = p_mid
            else:
                p_unstable = p_mid
                
        logger.info(f"[Dew] ✅ Bracket: unstable={p_unstable:.4f} -> stable={p_stable:.4f}")
        return p_unstable, p_stable

    def dp_process(self, lambd: float = 1.0) -> bool | None:
        cur_s_dp = self._find_S(self.p_i)
        while np.isclose(cur_s_dp['S'], 0.0, atol=1e-6):
            self.p_min_dew = self.p_i
            self.p_i = (self.p_max_dew + self.p_min_dew) / 2.0
            if (self.p_max_dew - self.p_min_dew) < 1e-12: return None
            cur_s_dp = self._find_S(self.p_i)

        s_dp_val = cur_s_dp['S']
        f_z, f_dp = cur_s_dp['letuch_z'], cur_s_dp['letuch_sp']
        r_dp = np.exp(f_z - f_dp) / s_dp_val
        y_dp = cur_s_dp['yi'] * np.power(r_dp, lambd)
        
        self._last_sum_y = float(np.sum(y_dp))
        ratio = np.clip(y_dp / self._z_arr, 1e-12, None)
        r_dp_safe = np.clip(r_dp, 1e-12, None)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio = np.log(ratio)
            safe_denom = np.where(np.abs(log_ratio) > 1e-12, log_ratio, 1.0)
            self.Sum = float(np.sum(np.log(r_dp_safe) / safe_denom))

        self._last_Ykz = float(np.sum(ratio))
        logger.debug(f"[Dew] P={self.p_i:.4f} | sum_y={self._last_sum_y:.6f} | Ykz²={self._last_Ykz**2:.2e}")

        if (abs(1.0 - self._last_sum_y) < 1e-4) or (self._last_Ykz ** 2 < 1e-4):
            return True
        else:
            self.p_max_dew = self.p_i  # 🔑 Сдвигаем верхнюю границу вниз (зеркально bubble point)
            self.p_i = (self.p_max_dew + self.p_min_dew) / 2.0
            return False

    def dp_convergence_loop(self, max_iter: int = 100) -> float | None:
        bracket = self._find_stability_bracket(self.p_min_dew, self.p_max_dew)
        if bracket[0] is not None:
            p_unstable, p_stable = bracket
            self.p_min_dew = p_unstable
            self.p_max_dew = p_stable
            self.p_i = (p_stable + p_unstable) / 2.0

        for i in range(1, max_iter + 1):
            if (self.p_max_dew - self.p_min_dew) < 1e-4:
                return self.p_i
            res = self.dp_process()
            if res is None or res:
                return self.p_i
        logger.warning(f"[Dew] ⚠️ Лимит итераций. P={self.p_i:.4f}")
        return self.p_i


# ==============================================================================
# 3. FAZOVAYA DIAGRAMMA (PHASE DIAGRAM)
# ==============================================================================
class PhaseDiagram:
    def __init__(self, composition: Composition, p_max_bar: float, t_min_c: float, t_max_c: float, t_step_c: float):
        self.composition = composition
        self.p_max = p_max_bar
        self.t_min = t_min_c + 273.15
        self.t_max = t_max_c + 273.15
        self.t_step = t_step_c
        self.temps = np.arange(self.t_min, self.t_max + self.t_step, self.t_step)
        self.results = {}

    def calc_phase_diagram(self, p_min_bar: float = 0.01):
        logger.info("🚀 Запуск расчёта фазовой диаграммы...")
        for t in self.temps:
            t_c = t - 273.15
            logger.info(f"🌡️ Обработка T={t_c:.2f}°C")
            
            # Bubble Point
            bp = SaturationPressure(self.composition, t, p_max_bar=self.p_max, p_min_bar=p_min_bar)
            pb = bp.sp_convergence_loop()
            
            # Dew Point
            dp = DewPointPressure(self.composition, t, p_max_bar=self.p_max, p_min_bar=p_min_bar)
            pdew = dp.dp_convergence_loop()
            
            self.results[t] = {'Pb': pb, 'Pdew': pdew}
            logger.info(f"✅ T={t_c:.2f}°C | Pb={pb:.4f} bar | Pdew={pdew:.4f} bar")
            
        return self.results

    def to_dataframe(self) -> pd.DataFrame:
        data = []
        for t, res in self.results.items():
            data.append({'T_C': t - 273.15, 'Pb_bar': res['Pb'], 'Pdew_bar': res['Pdew']})
        return pd.DataFrame(data)

    def plot(self, save_path: str | None = None):
        df = self.to_dataframe()
        plt.figure(figsize=(10, 6))
        
        plt.plot(df['T_C'], df['Pb_bar'], 'b-o', label='Bubble Point', markersize=4)
        plt.plot(df['T_C'], df['Pdew_bar'], 'r-o', label='Dew Point', markersize=4)
        
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Pressure (bar)')
        plt.title('Phase Envelope')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.xlim(df['T_C'].min() - 10, df['T_C'].max() + 10)
        
        if save_path:
            plt.savefig(save_path, dpi=300)
            logger.info(f"📊 График сохранён: {save_path}")
        plt.show()