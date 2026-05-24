import math
import logging
import numpy as np
from typing import Optional, Dict, Any

from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest
from _src.Composition.CompositionV2 import Composition

logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s | %(asctime)s | %(message)s', datefmt='%H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


class SaturationPressureCalculation:
    def __init__(self, composition_object: Composition, p_max: float, temp: float, p_min: float = 0.1):
        self.zi = composition_object
        # 🔑 Инициализация интервала. p_min не должен быть слишком маленьким (EOS ломается)
        self.p_min_bub = max(p_min, 0.5) 
        self.p_max_bub = p_max
        self.p_i = (self.p_max_bub + self.p_min_bub) / 2.0
        self.temp = temp + 273.14  # Точная константа из оригинала
        
        self.components = list(self.zi.composition.keys())
        self.z_dict = self.zi.composition

    def define_s_sp(self, p: float) -> Dict[str, Any]:
        """
        ТОЧНЫЙ ПОРТ ЛОГИКИ define_s_sp из VBA/эталона.
        Два независимых блока if, второй может перезаписать y_sp.
        """
        phase_stability = TwoPhaseStabilityTest(self.zi, p, self.temp)
        phase_stability.calculate_phase_stability()

        Sl = phase_stability.S_l
        Sv = phase_stability.S_v

        # V3: летучести исходной смеси (ln(f_i))
        letuch_z_arr = phase_stability._mixture_fugacities_arr
        letuch_z = {comp: float(letuch_z_arr[i]) for i, comp in enumerate(self.components)}

        # Тривиальное решение (однофазная система)
        if (Sl - 1.0) < 1e-5 and (Sv - 1.0) < 1e-5:
            y_sp = {comp: 0.0 for comp in self.components}
            return {'s_sp': 0.0, 'y_sp': y_sp, 'k_sp': None,
                    'letuch_sp': None, 'letuch_z': letuch_z}

        y_sp_dict = {}
        letuch_sp_dict = None
        k_sp = None

        # 🔑 БЛОК 1 (VBA: If Sl > 1 Then ... Else ... End If)
        if Sl > 1.0:
            if Sl > Sv:
                # Liquid branch: y = z / K_l
                k_sp = phase_stability.k_values_liquid
                letuch_sp_arr = phase_stability.liquid_eos.fugacities
                for comp in self.components:
                    y_sp_dict[comp] = self.z_dict[comp] / k_sp[comp]
            else:
                # Vapour branch: y = z * K_v
                k_sp = phase_stability.k_values_vapour
                letuch_sp_arr = phase_stability.vapour_eos.fugacities
                for comp in self.components:
                    y_sp_dict[comp] = self.z_dict[comp] * k_sp[comp]
            letuch_sp_dict = {comp: float(letuch_sp_arr[i]) for i, comp in enumerate(self.components)}
        else:
            # Sl <= 1
            if Sv < 1.0:
                y_sp_dict = {comp: 0.0 for comp in self.components}
                return {'s_sp': 0.0, 'y_sp': y_sp_dict, 'k_sp': None,
                        'letuch_sp': None, 'letuch_z': letuch_z}

        # 🔑 БЛОК 2 (VBA: If Sv > 1 Then ... Else ... End If)
        # Независимый блок, может перезаписать y_sp_dict и k_sp
        if Sv > 1.0:
            if Sv > Sl:
                # Vapour branch: y = z * K_v (перезапись!)
                k_sp = phase_stability.k_values_vapour
                letuch_sp_arr = phase_stability.vapour_eos.fugacities
                for comp in self.components:
                    y_sp_dict[comp] = self.z_dict[comp] * k_sp[comp]
                letuch_sp_dict = {comp: float(letuch_sp_arr[i]) for i, comp in enumerate(self.components)}
        else:
            # Sv <= 1
            if Sl < 1.0:
                y_sp_dict = {comp: 0.0 for comp in self.components}
                return {'s_sp': 0.0, 'y_sp': y_sp_dict, 'k_sp': None,
                        'letuch_sp': None, 'letuch_z': letuch_z}

        S_sp = sum(y_sp_dict.values())
        return {'s_sp': S_sp, 'y_sp': y_sp_dict, 'k_sp': k_sp,
                'letuch_sp': letuch_sp_dict, 'letuch_z': letuch_z}

    def calculate_saturation_pressure(self, lambd: float = 1.0, max_iter: int = 100) -> Optional[float]:
        """
        Устойчивый итерационный метод на основе эталонного кода.
        Использует брекетинг давления и термодинамические критерии сходимости.
        """
        tol_p = 1e-4      # Точность по давлению
        tol_s = 1e-4      # Точность по |1 - Sum(y)|
        tol_ykz_sq = 1e-4 # Точность по Ykz^2

        p_min = self.p_min_bub
        p_max = self.p_max_bub
        p_i = self.p_i

        # 🔑 ШАГ 0: Умная инициализация интервала
        # Проверяем, не находимся ли мы уже в тривиальной зоне, и расширяем интервал, если нужно
        res_min = self.define_s_sp(p_min)
        res_max = self.define_s_sp(p_max)
        
        # Если на обеих границах S_sp == 0 (тривиальное решение), интервал неверен
        if res_min['s_sp'] == 0.0 and res_max['s_sp'] == 0.0:
            logger.warning("Тривиальные решения на обеих границах. Корректировка интервала...")
            # Обычно P_sat находится там, где S_sp > 0. Попробуем уменьшить p_max
            while res_max['s_sp'] == 0.0 and p_max > p_min * 1.1:
                p_max /= 2.0
                res_max = self.define_s_sp(p_max)
            p_min = p_max / 10.0 # Сужаем нижнюю границу
            p_i = (p_min + p_max) / 2.0
            logger.info(f"Новый интервал: [{p_min:.3f}, {p_max:.3f}]")

        for iteration in range(max_iter):
            res = self.define_s_sp(p_i)
            s_sp = res['s_sp']

            # 🔑 Обработка тривиального решения (s_sp == 0)
            # Это значит, что давление слишком высокое (для bubble point) или низкое (для dew point)
            # Для bubble point: если S_sp=0, то мы выше точки насыщения -> уменьшаем p_max
            if s_sp == 0.0:
                p_max = p_i
                p_i = (p_max + p_min) / 2.0
                logger.debug(f"Iter {iteration:2d} | S_sp=0.0 | P_i={p_i:.4f} | Bracket=[{p_min:.3f}, {p_max:.3f}]")
                if p_max - p_min < tol_p:
                    logger.warning("Интервал сузился до tol_p при тривиальном решении.")
                    return p_i
                continue

            y_sp = res['y_sp']
            letuch_z = res['letuch_z']
            letuch_sp = res['letuch_sp']

            if letuch_sp is None:
                # Защита от некорректных данных
                p_min = p_i
                p_i = (p_max + p_min) / 2.0
                continue

            # 🔑 Расчёт r_sp (отношение летучестей)
            r_sp = {}
            for comp in self.components:
                # r_sp = f_z / (f_sp * S_sp)
                r_sp[comp] = math.exp(letuch_z[comp]) / (math.exp(letuch_sp[comp]) * s_sp)

            # 🔑 Обновление y_sp с релаксацией (lambda)
            y_sp_updated = {}
            for comp in self.components:
                y_sp_updated[comp] = y_sp[comp] * (r_sp[comp] ** lambd)

            sum_y_sp = sum(y_sp_updated.values())

            # 🔑 Расчёт Sum и Ykz для проверки сходимости
            Sum_val = 0.0
            Ykz = 0.0
            for comp in self.components:
                z_i = self.z_dict[comp]
                y_i = y_sp_updated[comp]
                r_i = r_sp[comp]
                
                # Защита от деления на ноль в логарифме
                ratio = y_i / z_i
                if ratio <= 0 or abs(ratio - 1.0) < 1e-12:
                    log_ratio = 1e-12
                else:
                    log_ratio = math.log(ratio)
                
                if abs(log_ratio) > 1e-15:
                    Sum_val += math.log(r_i) / log_ratio
                
                Ykz += y_i / z_i

            # 🔑 КРИТЕРИИ СХОДИМОСТИ (как в эталоне)
            # 1. |1 - Sum(y)| < tol_s
            # 2. Ykz^2 < tol_ykz_sq
            if abs(1.0 - sum_y_sp) < tol_s or (Ykz ** 2) < tol_ykz_sq:
                logger.info(f"✅ Сходимость на итерации {iteration}: P_sat = {p_i:.4f} bar | Σy={sum_y_sp:.5f} | Ykz²={Ykz**2:.2e}")
                return p_i

            # 🔑 Обновление брекета (бисекция, как в эталоне)
            # Если не сошлось, считаем, что мы ниже точки насыщения (для bubble point) -> увеличиваем p_min
            p_min = p_i
            p_i = (p_max + p_min) / 2.0
            
            logger.debug(f"Iter {iteration:2d} | Σy={sum_y_sp:.5f} | Ykz²={Ykz**2:.2e} | P_i={p_i:.4f} | Bracket=[{p_min:.3f}, {p_max:.3f}]")

            if p_max - p_min < tol_p:
                logger.info(f"Интервал сузился до {tol_p}. Возвращаем P_i={p_i:.4f}")
                return p_i

        logger.error(f"Превышено макс. число итераций ({max_iter}). Последний P_i={p_i:.4f}")
        return p_i

    def clear_cache(self):
        pass