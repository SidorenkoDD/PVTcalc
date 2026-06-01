from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton
import math as math
import numpy as np
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')


class SaturationPressure:
    def __init__(self, composition: Composition, t_K, p_max_bar:float = 1000, p_min_bar:float = 0.01):
        self.composition_object = composition
        self.t_K = t_K
        self.p_max_bar_sat = p_max_bar
        self.p_min_bar_sat = p_min_bar
        self.p_i = p_max_bar / 2
        self.p_max_bar_dew = p_max_bar
        self.p_min_bar_dew = p_min_bar

        self.zi = composition.composition
        self._components = tuple(self.zi.keys())
        self._components_mole_fractions = np.array(list(self.zi.values()))
        print(f'INIT_SUM_MOLE_FRACTIONS: {np.sum(self._components_mole_fractions)}')
        self._component_index = {comp: i for i, comp in enumerate(self._components)}
        self._nc = len(self._components)

    def _find_S(self, p):
        phase_stability_object = TwoPhaseStabilityTest(self.composition_object, p, t=self.t_K)
        phase_stability_object.calculate_phase_stability()
        print(f'S_v: {phase_stability_object.S_v}, S_l:{phase_stability_object.S_l}, stable: {phase_stability_object.stable}')
        ### Обе суммы < 1
        if ((phase_stability_object.S_l - 1) < (10 ** -5)) and ((phase_stability_object.S_v - 1) < (10 ** -5)):
            print('case0')
            yi = np.zeros(self._nc)
            return {'S': 0, 'yi': yi, 'k_sp': None,
                        'letuch_sp': None, 'letuch_z': None}
        
        ### Обе фазы больше 1: двухфазная область
        if phase_stability_object.S_l > 1:
            ### жидкая фаза менее стабильна
            if phase_stability_object.S_l > phase_stability_object.S_v:
                print('case1')
                k_values = phase_stability_object._k_l_arr.copy()
                fugacities_mix = phase_stability_object._mixture_fugacities_arr
                fugacities_phase = phase_stability_object.liquid_eos.fugacities
                yi = self._components_mole_fractions / k_values
            else:
                ### газовая фаза менее стабильна
                print('case2')
                k_values = phase_stability_object._k_v_arr.copy()
                fugacities_mix = phase_stability_object._mixture_fugacities_arr
                fugacities_phase = phase_stability_object.vapour_eos.fugacities
                yi = self._components_mole_fractions * k_values

        else:
            if phase_stability_object.S_v < 1:
                print('case3')
                yi = np.zeros(self._nc)
                return {'S': 0, 'yi': yi, 'k_sp': None,
                        'letuch_sp': None, 'letuch_z': None}

        if phase_stability_object.S_v > 1:
            if phase_stability_object.S_v > phase_stability_object.S_l:
                print('case4')
                k_values = phase_stability_object._k_v_arr.copy()
                fugacities_mix = phase_stability_object._mixture_fugacities_arr
                fugacities_phase = phase_stability_object.vapour_eos.fugacities
                yi = self._components_mole_fractions * k_values

            else:
                if phase_stability_object.S_l < 1:
                    print('case5')
                    yi = np.zeros(self._nc)
                    return {'S': 0, 'yi': yi, 'k_sp': None,
                            'letuch_sp': None, 'letuch_z': None}

        sum_mole_fractions = sum(yi)
        return {'S': sum_mole_fractions, 'yi': yi, 'k_sp': k_values,  
                'letuch_sp': fugacities_phase, 'letuch_z': fugacities_mix}

    def sp_process(self, lambd=1):
        cur_s_sp = self._find_S(self.p_i)
        print(f'S: {cur_s_sp['S']}')
        # 1. Безопасная проверка на ноль
        while np.isclose(cur_s_sp['S'], 0.0, atol= 1e-6):
            self.p_max_bar_sat = self.p_i
            self.p_i = (self.p_max_bar_sat + self.p_min_bar_sat) / 2.0

            if (self.p_max_bar_sat - self.p_min_bar_sat) < 1e-12:
                print('Строка 70: sp_process не сошелся. В результате None')
                return None
                
            cur_s_sp = self._find_S(self.p_i)
            
        s_sp_val = cur_s_sp['S'] 
        
        # 3. Все math.* заменяем на np.* (векторные операции)
        # exp_z = np.exp(cur_s_sp['letuch_z'])
        # exp_sp = np.exp(cur_s_sp['letuch_sp'])
        
        # r_sp = exp_z / (exp_sp * s_sp_val)  
        log_ratio_fug = cur_s_sp['letuch_z'] - cur_s_sp['letuch_sp']
        r_sp = np.exp(log_ratio_fug) / s_sp_val

        y_sp = cur_s_sp['yi'] * np.power(r_sp, lambd)

        self.sum_y_sp = np.sum(y_sp)
        print(f'S: {self.sum_y_sp}')
        # 4. Логарифмы и суммы поэлементно
        ratio = y_sp / self._components_mole_fractions

        # Защитное ограничение снизу (избегаем log(0) и деления на 0)
        ratio = np.clip(ratio, 1e-12, None)
        r_sp = np.clip(r_sp, 1e-12, None)

        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio = np.log(ratio)
            # Заменяем 0 в знаменателе на nan, чтобы не ломать сумму, или используем safe_div
            safe_div = np.where(np.abs(log_ratio) > 1e-12, np.log(r_sp) / log_ratio, 0.0)
            self.Sum = np.sum(safe_div)

        self.Ykz = np.sum(ratio)

        print(f"p={self.p_i:.2f}, sum_y_sp={self.sum_y_sp:.6f}, Ykz={self.Ykz:.6f}")
        # 5. Проверка сходимости
        if (np.abs(1.0 - self.sum_y_sp) < 1e-4) or (self.Ykz ** 2 < 1e-4):
            print(f'Pb найдено: {self.p_i}')
        else:
            self.p_min_bar_sat = self.p_i
            self.p_i = (self.p_max_bar_sat + self.p_min_bar_sat) / 2.0

    def sp_convergence_loop(self):
        self.sp_process()
        if self.p_max_bar_sat - self.p_min_bar_sat < math.pow(10, -4):
            return None
        
        while not (abs(1 - self.sum_y_sp) < math.pow(10, -4) or math.pow(self.Ykz, 2) < math.pow(10, -4)):
            self.sp_process()
            if self.p_max_bar_sat - self.p_min_bar_sat < math.pow(10, -4):
                return None

        self.p_b = self.p_i
        self.p_i = self.p_i / 2
        return self.p_b
    

class SaturationPressureFromFlash:
    def __init__(self, composition: Composition, t_K: float, p_max_bar: float = 1000, p_min_bar:float = 0.1):
        self.composition = composition
        self.t_K = t_K
        self.p_max_bar = p_max_bar
        self.p_min_bar = p_min_bar
        self.p_i = self.p_max_bar/2

    def loop(self):
        phase_stability_object = TwoPhaseStabilityTest(composition=self.composition, p= self.p_i, t = self.t_K)
        phase_stability_object.calculate_phase_stability()

        while (phase_stability_object.stable) == True:
            self.p_max_bar = self.p_i
            self.p_i = (self.p_max_bar+self.p_min_bar) / 2
            print(self.p_i)
            phase_stability_object = TwoPhaseStabilityTest(composition=self.composition, p= self.p_i, t = self.t_K)
            phase_stability_object.calculate_phase_stability()

        phase_equil_object = PhaseEquilibriumNewton(self.composition, p = self.p_i, t = self.t_K, k_values=phase_stability_object.k_values_for_flash)
        flash_result = phase_equil_object.find_solve_loop()
        print(flash_result['Fl'])
        print(f'Pi: {self.p_i}')
        print(f'Pmax {self.p_max_bar}')
        print(f'Pmin {self.p_min_bar}')

        while (flash_result['Fl'] < 0.99):
            self.p_min_bar = self.p_i
            self.p_i = (self.p_max_bar + self.p_min_bar) / 2
            print(f'Upd p_i: {self.p_i}')
            phase_stability_object = TwoPhaseStabilityTest(composition=self.composition, p= self.p_i, t = self.t_K)
            phase_stability_object.calculate_phase_stability()
            print(f'phase_stab_upd_pi: {phase_stability_object.stable}')
            if phase_stability_object.stable == True:
                self.p_max_bar = self.p_i
                self.p_i = (self.p_max_bar + self.p_min_bar) / 2
            else:
                phase_equil_object = PhaseEquilibriumNewton(self.composition, p = self.p_i, t = self.t_K, k_values=phase_stability_object.k_values_for_flash)
                flash_result = phase_equil_object.find_solve_loop()
                print(flash_result['Fl'])
                print(f'Pi: {self.p_i}')
                print(f'Pmax {self.p_max_bar}')
                print(f'Pmin {self.p_min_bar}')
                print('======')

        self.p_bub = self.p_i

        print(self.p_bub)


