import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest


class SaturationPressure:
    def __init__(self, composition: Composition, p_max: float, temp: float, p_min: float = 0.1):
        self.composition = composition
        self.zi = composition.composition  # dict мольных долей
        self.p_min_bub = p_min
        self.p_max_bub = p_max
        self.p_i = self.p_max_bub / 2
        self.p_min_dew = 0.0
        self.p_max_dew = p_max
        self.temp = temp

    def define_s_sp(self, p: float):
        phase_stability = TwoPhaseStabilityTest(composition=self.composition, p=p, t=self.temp)
        phase_stability.calculate_phase_stability()
        print(f'phase stab Sv: {phase_stability.S_v}, Sl: {phase_stability.S_l}, stable: {phase_stability.stable}')
        # Летучести исходной смеси (аналог phase_stability.initial_eos.fugacities)
        comps = list(self.zi.keys())
        letuch_z = {c: float(phase_stability._mixture_fugacities_arr[i]) for i, c in enumerate(comps)}

        # Проверка на стабильность / тривиальное решение
        if (phase_stability.S_l - 1) < 1e-5 and (phase_stability.S_v - 1) < 1e-5:
            y_sp = {c: 0 for c in comps}
            return {'s_sp': 0, 'y_sp': y_sp, 'k_sp': None, 'r_sp': None, 
                    'letuch_sp': None, 'letuch_z': letuch_z}

        letuch_l = {c: float(phase_stability.liquid_eos.fugacities[i]) for i, c in enumerate(comps)}
        letuch_v = {c: float(phase_stability.vapour_eos.fugacities[i]) for i, c in enumerate(comps)}

        # --- Блок 1 (1:1 с PhaseDiagram_v4.py) ---
        if phase_stability.S_l > 1:
            if phase_stability.S_l > phase_stability.S_v:
                k_sp = phase_stability.k_values_liquid
                r_sp = phase_stability.ri_l  # ri_l не экспортируется в V3, но в sp_process он всё равно пересчитывается заново
                letuch_sp = letuch_l
                y_sp = {c: self.zi[c] / phase_stability.k_values_liquid[c] for c in comps}
            else:
                k_sp = phase_stability.k_values_vapour
                r_sp = phase_stability.ri_v
                letuch_sp = letuch_v
                y_sp = {c: self.zi[c] * phase_stability.k_values_vapour[c] for c in comps}
        else:
            if phase_stability.S_v < 1:
                y_sp = {c: 0 for c in comps}
                return {'s_sp': 0, 'y_sp': y_sp, 'k_sp': None, 'r_sp': None, 
                        'letuch_sp': None, 'letuch_z': letuch_z}

        # --- Блок 2 (1:1 с PhaseDiagram_v4.py) ---
        if phase_stability.S_v > 1:
            if phase_stability.S_v > phase_stability.S_l:
                k_sp = phase_stability.k_values_vapour
                r_sp = phase_stability.ri_v
                letuch_sp = letuch_v
                y_sp = {c: self.zi[c] * phase_stability.k_values_vapour[c] for c in comps}
            else:
                if phase_stability.S_l < 1:
                    y_sp = {c: 0 for c in comps}
                    return {'s_sp': 0, 'y_sp': y_sp, 'k_sp': None, 'r_sp': None, 
                            'letuch_sp': None, 'letuch_z': letuch_z}

        S_sp = sum(y_sp.values())
        return {'s_sp': S_sp, 'y_sp': y_sp, 'k_sp': k_sp, 'r_sp': r_sp, 
                'letuch_sp': letuch_sp, 'letuch_z': letuch_z}

    def sp_process(self, lambd=1):
        cur_s_sp = self.define_s_sp(self.p_i)
        print(cur_s_sp)
        # Если s_sp == 0, сужаем диапазон в сторону уменьшения давления
        while cur_s_sp['s_sp'] == 0:
            self.p_max_bub = self.p_i
            self.p_i = (self.p_max_bub + self.p_min_bub) / 2
            if self.p_max_bub - self.p_min_bub < 1e-5:
                return None
            cur_s_sp = self.define_s_sp(self.p_i)

        # Пересчёт r_sp и y_sp (строго как в оригинале)
        r_sp = {}
        for component in cur_s_sp['letuch_z'].keys():
            r_sp[component] = math.exp(cur_s_sp['letuch_z'][component]) / (
                math.exp(cur_s_sp['letuch_sp'][component]) * cur_s_sp['s_sp']
            )

        y_sp = {component: cur_s_sp['y_sp'][component] * math.pow(r_sp[component], lambd) 
                for component in r_sp.keys()}

        self.sum_y_sp = sum(y_sp.values())
        self.Sum = sum(math.log(r_sp[component]) / math.log(y_sp[component] / self.zi[component]) 
                       for component in self.zi.keys())
        self.Ykz = sum(y_sp[component] / self.zi[component] for component in self.zi.keys())

        # Проверка сходимости
        if (abs(1 - self.sum_y_sp) < 1e-3) or (math.pow(self.Ykz, 2) < 1e-4):
            print(f'Pb найдено: {self.p_i}')
        else:
            self.p_min_bub = self.p_i
            self.p_i = (self.p_max_bub + self.p_min_bub) / 2

    def sp_convergence_loop(self):
        self.sp_process()
        if self.p_max_bub - self.p_min_bub < 1e-5:
            return None
        while not (abs(1 - self.sum_y_sp) < 1e-4 or math.pow(self.Ykz, 2) < 1e-4):
            self.sp_process()
            if self.p_max_bub - self.p_min_bub < 1e-5:
                return None
        self.p_b = self.p_i
        self.p_i = self.p_i / 2  # Подготовка к dew-point
        return self.p_b

    # ==================== DEW POINT (аналогично Bubble Point) ====================

    def define_s_dp(self, p: float):
        phase_stability = TwoPhaseStabilityTest(composition=self.composition, p=p, t=self.temp)
        phase_stability.calculate_phase_stability()

        comps = list(self.zi.keys())
        letuch_z = {c: float(phase_stability._mixture_fugacities_arr[i]) for i, c in enumerate(comps)}

        if (phase_stability.S_l - 1) < 1e-5 and (phase_stability.S_v - 1) < 1e-5:
            y_dp = {c: 0 for c in comps}
            return {'s_dp': 0, 'y_dp': y_dp, 'k_dp': None, 'r_dp': None, 
                    'letuch_dp': None, 'letuch_z': letuch_z}

        letuch_l = {c: float(phase_stability.liquid_eos.fugacities[i]) for i, c in enumerate(comps)}
        letuch_v = {c: float(phase_stability.vapour_eos.fugacities[i]) for i, c in enumerate(comps)}

        if phase_stability.S_l > 1:
            if phase_stability.S_l > phase_stability.S_v:
                k_dp = phase_stability.k_values_liquid
                r_dp = None
                letuch_dp = letuch_l
                y_dp = {c: self.zi[c] / phase_stability.k_values_liquid[c] for c in comps}
            else:
                k_dp = phase_stability.k_values_vapour
                r_dp = None
                letuch_dp = letuch_v
                y_dp = {c: self.zi[c] * phase_stability.k_values_vapour[c] for c in comps}
        else:
            if phase_stability.S_v < 1:
                y_dp = {c: 0 for c in comps}
                return {'s_dp': 0, 'y_dp': y_dp, 'k_dp': None, 'r_dp': None, 
                        'letuch_dp': None, 'letuch_z': letuch_z}

        if phase_stability.S_v > 1:
            if phase_stability.S_v > phase_stability.S_l:
                k_dp = phase_stability.k_values_vapour
                r_dp = None
                letuch_dp = letuch_v
                y_dp = {c: self.zi[c] * phase_stability.k_values_vapour[c] for c in comps}
            else:
                if phase_stability.S_l < 1:
                    y_dp = {c: 0 for c in comps}
                    return {'s_dp': 0, 'y_dp': y_dp, 'k_dp': None, 'r_dp': None, 
                            'letuch_dp': None, 'letuch_z': letuch_z}

        S_dp = sum(y_dp.values())
        return {'s_dp': S_dp, 'y_dp': y_dp, 'k_dp': k_dp, 'r_dp': r_dp, 
                'letuch_dp': letuch_dp, 'letuch_z': letuch_z}

    def dp_process(self, lambd=1):
        cur_s_dp = self.define_s_dp(self.p_i)

        while cur_s_dp['s_dp'] == 0:
            self.p_min_dew = self.p_i
            self.p_i = (self.p_max_dew + self.p_min_dew) / 2
            if self.p_max_dew - self.p_min_dew < 1e-5:
                return None
            cur_s_dp = self.define_s_dp(self.p_i)

        r_dp = {}
        for component in cur_s_dp['letuch_z'].keys():
            r_dp[component] = math.exp(cur_s_dp['letuch_z'][component]) / (
                math.exp(cur_s_dp['letuch_dp'][component]) * cur_s_dp['s_dp']
            )

        y_dp = {component: cur_s_dp['y_dp'][component] * math.pow(r_dp[component], lambd) 
                for component in r_dp.keys()}

        self.sum_y_dp = sum(y_dp.values())
        self.Sum_dp = sum(math.log(r_dp[component]) / math.log(y_dp[component] / self.zi[component]) 
                          for component in self.zi.keys())
        self.Ykz_dp = sum(y_dp[component] / self.zi[component] for component in self.zi.keys())

        if abs(1 - self.sum_y_dp) < 1e-3 or math.pow(self.Ykz_dp, 2) < 1e-3:
            print(f'Pdew найдено: {self.p_i}')
        else:
            self.p_max_dew = self.p_i
            self.p_i = (self.p_max_dew + self.p_min_dew) / 2

    def dp_convergence_loop(self):
        self.dp_process()
        if self.p_max_dew - self.p_min_dew < 1e-5:
            return None
        while not (abs(1 - self.sum_y_dp) < 1e-3 or math.pow(self.Ykz_dp, 2) < 1e-3):
            self.dp_process()
            if self.p_max_dew - self.p_min_dew < 1e-5:
                return None
        self.p_dew = self.p_i
        return self.p_dew


class PhaseDiagram:
    def __init__(self, composition: Composition, p_max: float, t_min: float, t_max: float, t_step: float, c7_plus_correlations: dict = None):
        self.composition = composition
        # Свойства компонент вычисляются ОДИН раз до цикла

        
        self.p_max = p_max
        self.t_min = t_min
        self.t_max = t_max
        self.t_step = t_step
        self.temp_arange = np.arange(self.t_min + 273.15, self.t_max + 273.15, self.t_step)
        self.results = {}

    def calc_phase_diagram(self):
        for temp in self.temp_arange:
            cur_sat = SaturationPressure(self.composition, self.p_max, temp)
            pb = cur_sat.sp_convergence_loop()
            pdew = cur_sat.dp_convergence_loop()
            self.results[temp] = [pb, pdew]
            print(f"Temp: {temp-273.15:.2f}°C, Pb: {pb}, Pdew: {pdew}")

    def plot_phase_diagram(self):
        bubble_points = []
        dew_points = []
        temps = []
        
        for temp, (pb, pdew) in self.results.items():
            temps.append(temp - 273.15)
            bubble_points.append(pb if pb is not None else np.nan)
            dew_points.append(pdew if pdew is not None else np.nan)
        
        plt.figure(figsize=(10, 6))
        plt.plot(temps, np.array(bubble_points), 'b-o', label='Bubble Point')
        plt.plot(temps, np.array(dew_points), 'r-o', label='Dew Point')
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Pressure (MPa)')
        plt.title('Phase Diagram')
        plt.legend()
        plt.grid(True)
        plt.ylim(0, self.p_max)
        plt.show()

    def get_phase_diagram_data(self):
        temps, bps, dps = [], [], []
        for temp, (pb, pdew) in self.results.items():
            temps.append(temp - 273.15)
            bps.append(pb if pb is not None else np.nan)
            dps.append(pdew if pdew is not None else np.nan)
        return pd.DataFrame({'Temp': temps, 'Bubble': bps, 'Dew': dps})