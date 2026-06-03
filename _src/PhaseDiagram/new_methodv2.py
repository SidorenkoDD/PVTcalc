import numpy as np
import logging
from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(name)s | %(levelname)s | %(message)s')

class SaturationPressure:
    def __init__(self, composition: Composition, t_K: float, p_max_bar: float = 1000, p_min_bar: float = 0.01):
        self.composition = composition
        self.composition.T = t_K
        #self.composition_object = composition.new_composition(composition.composition, deep_copy=True)
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
        logger.info(f"Инициализация: T={self.t_K:.2f} K, P_min={self.p_min_bar_sat:.4f}, P_max={self.p_max_bar_sat:.4f}")

    def _find_S(self, p: float) -> dict:
        """Точная копия логики выбора фазы из PhaseDiagram_v4.py, векторизованная."""
        phase_stability = TwoPhaseStabilityTest(self.composition, p, t=self.t_K)
        phase_stability.calculate_phase_stability()
        
        S_l, S_v = phase_stability.S_l, phase_stability.S_v
        
        # Case 0: Однофазная стабильная система
        if (S_l - 1.0) < 1e-5 and (S_v - 1.0) < 1e-5:
            return {'S': 0.0, 'yi': np.zeros(self._nc), 'letuch_sp': None, 'letuch_z': None, 'stable': phase_stability.stable}

        yi = np.zeros(self._nc)
        letuch_z = phase_stability._mixture_fugacities_arr.copy()
        letuch_sp = None

        # Строгая последовательность условий как в VBA/v4
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

        S_sp = float(np.sum(yi))
        return {'S': S_sp, 'yi': yi, 'letuch_sp': letuch_sp, 'letuch_z': letuch_z, 'stable': phase_stability.stable}
    

    def _find_stability_bracket(self, p_low_init=0.01, p_high_init=1000.0) -> tuple[float, float]:
        """
        Исправленная стратегия:
        1. Проверяем статус на границах.
        2. Если обе в одном состоянии — ищем переход в нужном направлении.
        3. Сужаем бисекцией ТОЛЬКО если нашли разные статусы.
        """
        logger.info(f"🔍 Поиск bracket: start=[{p_low_init:.3f}, {p_high_init:.3f}]")

        # 1. Проверяем начальные точки
        res_low = self._find_S(p_low_init)
        res_high = self._find_S(p_high_init)
        
        stable_low = res_low['stable']
        stable_high = res_high['stable']
        
        # 2. Если статусы одинаковые — ищем переход
        if stable_low == stable_high:
            if stable_high:  # Обе стабильны → ищем unstable ниже
                p_search = p_high_init
                for step in range(10):
                    p_search *= 0.5
                    if p_search < p_low_init * 0.1:
                        logger.warning("Не найдена unstable точка при снижении")
                        return p_low_init, p_high_init
                    if not self._find_S(p_search)['stable']:
                        p_unstable = p_search
                        p_stable = p_search * 2.0  # Возвращаемся на шаг вверх
                        break
                else:
                    logger.warning("Превышен лимит поиска unstable")
                    return p_low_init, p_high_init
            else:  # Обе нестабильны → ищем stable выше
                p_search = p_high_init
                for step in range(10):
                    p_search *= 1.5
                    if p_search > 3000:
                        logger.warning("Не найдена stable точка при повышении")
                        return p_low_init, p_high_init
                    if self._find_S(p_search)['stable']:
                        p_stable = p_search
                        p_unstable = p_search / 1.5  # Возвращаемся на шаг вниз
                        break
                else:
                    logger.warning("Превышен лимит поиска stable")
                    return p_low_init, p_high_init
        else:
            # Статусы разные — определяем, какая граница какая
            if stable_low:
                p_stable, p_unstable = p_low_init, p_high_init
            else:
                p_stable, p_unstable = p_high_init, p_low_init

        # 3. Гарантируем порядок: stable > unstable для bubble point
        if p_stable < p_unstable:
            p_stable, p_unstable = p_unstable, p_stable

        # 4. Сужаем интервал бисекцией (но с проверкой, что он не нулевой)
        if abs(p_stable - p_unstable) < 1e-6:
            logger.warning(f"Интервал слишком мал после поиска: {p_stable:.4f} vs {p_unstable:.4f}")
            # Расширяем искусственно для дальнейшего поиска
            p_unstable *= 0.99
            p_stable *= 1.01

        for i in range(10):  # 10 шагов = сужение в 1024 раза
            p_mid = (p_stable + p_unstable) / 2.0
            if self._find_S(p_mid)['stable']:
                p_stable = p_mid
            else:
                p_unstable = p_mid
            # Защита от схлопывания
            if abs(p_stable - p_unstable) < 1e-8:
                break
                
        logger.info(f"Диапазон поиска: unstable={p_unstable:.4f} bar -> stable={p_stable:.4f} bar")
        return p_unstable, p_stable

    def sp_process(self, lambd: float = 1.0) -> bool | None:
        cur_s_sp = self._find_S(self.p_i)

        # Если система однофазная, сужаем p_max вниз (как в оригинале)
        while np.isclose(cur_s_sp['S'], 0.0, atol=1e-6):
            self.p_max_bar_sat = self.p_i
            self.p_i = (self.p_max_bar_sat + self.p_min_bar_sat) / 2.0
            if (self.p_max_bar_sat - self.p_min_bar_sat) < 1e-12:
                return None
            cur_s_sp = self._find_S(self.p_i)

        s_sp_val = cur_s_sp['S']
        f_z = cur_s_sp['letuch_z']
        f_sp = cur_s_sp['letuch_sp']

        # Формула Rsp из оригинала
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
        logger.debug(f"P={self.p_i:.4f} | sum_y={self._last_sum_y:.6f} | Ykz²={self._last_Ykz**2:.2e}")

        # Критерий сходимости
        if (abs(1.0 - self._last_sum_y) < 1e-4) or (self._last_Ykz ** 2 < 1e-4):
            return True
        else:
            # 🔑 ИСПРАВЛЕНИЕ: двусторонняя бисекция вместо односторонней
            if self._last_sum_y > 1.0:
                # sum_y > 1 → мы в двухфазной области → нужно повысить давление
                self.p_min_bar_sat = self.p_i
            else:
                # sum_y < 1 → мы в однофазной области → нужно понизить давление
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
            logger.info(f"Старт бисекции с P={self.p_i:.4f} bar, интервал={p_stable - p_unstable:.4f} bar")
        else:
            logger.warning("Границы не найдены или нулевые, используем исходные границы")
            # Если bracket не удался, пробуем хотя бы запустить с исходными границами
            self.p_min_bar_sat = 0.01
            self.p_max_bar_sat = 1000.0
            self.p_i = 500.0
            
        # Запускаем цикл с защитой от зацикливания
        prev_p = None
        for i in range(1, max_iter + 1):
            # Проверка на "залипание"
            if prev_p is not None and abs(self.p_i - prev_p) < 1e-10:
                logger.warning(f"Давление застряло на итерации {i}. P={self.p_i:.6f}")
                break
            prev_p = self.p_i
            
            if (self.p_max_bar_sat - self.p_min_bar_sat) < 1e-6:
                logger.info(f"Интервал схлопнулся на итерации {i}. Возврат P={self.p_i:.4f} bar")
                return self.p_i
                
            res = self.sp_process()
            if res is None:
                logger.warning(f"sp_process вернул None на итерации {i}")
                return self.p_i
            if res:
                logger.info(f"Сходимость достигнута на итерации {i}. P_sat={self.p_i:.4f} bar")
                return self.p_i
                
        logger.warning(f" Превышен лимит итераций ({max_iter}). Последнее P={self.p_i:.4f} bar, интервал={self.p_max_bar_sat - self.p_min_bar_sat:.4e}")
        return self.p_i

class SaturationPressureFromFlash:
    def __init__(self, composition: Composition, t_K: float, p_max_bar: float = 1000, p_min_bar:float = 0.1):
        self.composition = composition
        self.t_K = t_K
        self.p_max_bar = p_max_bar
        self.p_min_bar = p_min_bar
        self.p_i = self.p_max_bar / 2.0

    def loop(self):
        phase_stability_object = TwoPhaseStabilityTest(composition=self.composition, p=self.p_i, t=self.t_K)
        phase_stability_object.calculate_phase_stability()

        while phase_stability_object.stable:
            self.p_max_bar = self.p_i
            self.p_i = (self.p_max_bar + self.p_min_bar) / 2.0
            logger.debug(f"Flash stable=True, new P={self.p_i:.4f}")
            phase_stability_object = TwoPhaseStabilityTest(composition=self.composition, p=self.p_i, t=self.t_K)
            phase_stability_object.calculate_phase_stability()

        from _src.VLE.PhaseEquilibriumNewtonV2 import PhaseEquilibriumNewton
        phase_equil_object = PhaseEquilibriumNewton(self.composition, p=self.p_i, t=self.t_K, k_values=phase_stability_object.k_values_for_flash)
        flash_result = phase_equil_object.find_solve_loop()
        logger.debug(f"Flash Fl={flash_result['Fl']:.6f} at P={self.p_i:.4f}")

        while flash_result['Fl'] < 0.99:
            self.p_min_bar = self.p_i
            self.p_i = (self.p_max_bar + self.p_min_bar) / 2.0
            phase_stability_object = TwoPhaseStabilityTest(composition=self.composition, p=self.p_i, t=self.t_K)
            phase_stability_object.calculate_phase_stability()
            
            if phase_stability_object.stable:
                self.p_max_bar = self.p_i
                self.p_i = (self.p_max_bar + self.p_min_bar) / 2.0
            else:
                phase_equil_object = PhaseEquilibriumNewton(self.composition, p=self.p_i, t=self.t_K, k_values=phase_stability_object.k_values_for_flash)
                flash_result = phase_equil_object.find_solve_loop()
                logger.debug(f"Upd P={self.p_i:.4f}, Fl={flash_result['Fl']:.6f}")

        self.p_bub = self.p_i
        logger.info(f"✅ Bubble Point from Flash найден: P={self.p_bub:.4f} bar")
        return self.p_bub