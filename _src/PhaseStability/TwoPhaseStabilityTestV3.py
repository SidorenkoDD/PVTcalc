import numpy as np
from _src.EOS.BrusilovskiyEOSV2 import BrusilovskiyEOS
from _src.Composition.CompositionV2 import Composition
from _src.Utils.Errors import StopIterationError
from _src.Utils.Constants import (
    TOL_TWO_PHASE_STABILITY_CONVERGENCE,
    TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION,
)
import logging

logger = logging.getLogger('MBALPVT.PVTDataModel.PVTCore.TwoPhaseStabilityTest')


class TwoPhaseStabilityTest:
    """
    Тест фазовой стабильности по методу Михельсена.
    Определяет, является ли смесь однофазной или может распасться на две фазы.
    """
    
    def __init__(self, composition: Composition, p: float, t: float):
        self.p = float(p)
        self.t = float(t)
        
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
        
        # Векторное представление
        self._components = tuple(self.zi.keys())
        self._component_index = {comp: i for i, comp in enumerate(self._components)}
        self._nc = len(self._components)
        
        self._z_feed = np.fromiter(
            (self.zi[c] for c in self._components), dtype=np.float64, count=self._nc
        )
        self._acentric_factor = np.fromiter(
            (self.composition_data['acentric_factor'][c] for c in self._components),
            dtype=np.float64, count=self._nc
        )
        self._critical_pressure = np.fromiter(
            (self.composition_data['critical_pressure'][c] for c in self._components),
            dtype=np.float64, count=self._nc
        )
        self._critical_temperature = np.fromiter(
            (self.composition_data['critical_temperature'][c] for c in self._components),
            dtype=np.float64, count=self._nc
        )
        
        self._k_v_arr = None
        self._k_l_arr = None
        self._k_flash_arr = None
        self._yi_v_arr = None
        self._xi_l_arr = None
        self._mixture_fugacities_arr = None
        
        # Публичные dict-представления
        self.yi_v = None
        self.xi_l = None

    # =====================================================================
    # ВСПОМОГАТЕЛЬНЫЕ МЕТОДЫ
    # =====================================================================
    
    def _array_to_dict(self, arr: np.ndarray) -> dict:
        return {comp: float(arr[i]) for i, comp in enumerate(self._components)}
    
    def _dict_to_array(self, data: dict) -> np.ndarray:
        return np.fromiter((data[c] for c in self._components), dtype=np.float64, count=self._nc)
    
    @staticmethod
    def _safe_log_metric(arr: np.ndarray) -> float:
        return float(np.sum(np.log(arr) ** 2))

    # =====================================================================
    # WILSON K-VALUES
    # =====================================================================
    
    def _calc_k_values_wilson_array(self) -> np.ndarray:
        return np.exp(
            5.37 * (1.0 + self._acentric_factor) * (1.0 - (self._critical_temperature / self.t))
        ) / (self.p / self._critical_pressure)

    # =====================================================================
    # РАСЧЕТ ПРОБНЫХ СОСТАВОВ
    # =====================================================================
    
    def _calc_Yi_v(self, k_values_vapour: np.ndarray) -> np.ndarray:
        return self._z_feed * k_values_vapour
    
    def _calc_Xi_l(self, k_values_liquid: np.ndarray) -> np.ndarray:
        return self._z_feed / k_values_liquid
    
    @staticmethod
    def _calc_S(arr: np.ndarray) -> float:
        return float(np.sum(arr))
    
    @staticmethod
    def _normalize_mole_fractions(arr: np.ndarray, S: float) -> np.ndarray:
        return arr / S if S > 1e-15 else arr

    # =====================================================================
    # Ri - обновление K-факторов
    # =====================================================================
    
    @staticmethod
    def _calc_ri_vapour(vapour_fugacities: np.ndarray, mixture_fugacities: np.ndarray, S_v: float) -> np.ndarray:
        return np.exp(mixture_fugacities - vapour_fugacities) / S_v
    
    @staticmethod
    def _calc_ri_liquid(liquid_fugacities: np.ndarray, mixture_fugacities: np.ndarray, S_l: float) -> np.ndarray:
        return np.exp(liquid_fugacities - mixture_fugacities) * S_l
    
    @staticmethod
    def _update_k_values(k_values_old: np.ndarray, ri: np.ndarray) -> np.ndarray:
        return k_values_old * ri
    
    @staticmethod
    def _check_convergence(ri: np.ndarray) -> bool:
        return float(np.sum((ri - 1.0) ** 2)) < TOL_TWO_PHASE_STABILITY_CONVERGENCE
    
    @staticmethod
    def _check_trivial_solution(k_values: np.ndarray) -> bool:
        return float(np.sum(np.log(k_values) ** 2)) < TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION

    # =====================================================================
    # ИНТЕРПРЕТАЦИЯ РЕЗУЛЬТАТОВ
    # =====================================================================
    
    def _interpretate_stability_analysis(self):
        """
        Определяет тип фазового поведения на основе результатов двух тестов.
        """
        # Случай 1: Смесь стабильна (однофазная)
        if (self.stable is None) or (
            (self.S_v <= 1.0) and (self.S_l <= 1.0) or
            ((self.S_v <= 1.0) and self.convergence_trivial_solution_l) or
            (self.convergence_trivial_solution_v and (self.S_l <= 1.0)) or
            (self.convergence_trivial_solution_v and self.convergence_trivial_solution_l)
        ):
            self.stable = True
            self._k_flash_arr = None
            self.k_values_for_flash = None
            logger.debug(f"Стабильная смесь: S_v={self.S_v:.3f}, S_l={self.S_l:.3f}")
            return
        
        # Случай 2: Только vapour-тест показывает нестабильность → BUBBLE POINT
        if (self.S_v > 1.0) and (self.S_l <= 1.0):
            self.stable = False
            self._k_flash_arr = self._k_v_arr.copy()
            self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)
            logger.debug(f"Bubble point: S_v={self.S_v:.3f} > 1, S_l={self.S_l:.3f} <= 1")
            return
        
        # Случай 3: Только liquid-тест показывает нестабильность → DEW POINT
        if (self.S_l > 1.0) and (self.S_v <= 1.0):
            self.stable = False
            self._k_flash_arr = self._k_l_arr.copy()
            self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)
            logger.debug(f"Dew point: S_l={self.S_l:.3f} > 1, S_v={self.S_v:.3f} <= 1")
            return
        
        # Случай 4: Оба теста показывают нестабильность → выбираем по соотношению
        if (self.S_v > 1.0) and (self.S_l > 1.0):
            ratio = self.S_v / self.S_l if self.S_l > 1e-10 else float('inf')
            
            if ratio > 1.5:
                # vapour-тест доминирует → ближе к bubble point
                self.stable = False
                self._k_flash_arr = self._k_v_arr.copy()
                self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)
                logger.debug(f"Двухфазная (bubble-наклон): S_v/S_l={ratio:.2f}")
            elif ratio < 0.67:
                # liquid-тест доминирует → ближе к dew point
                self.stable = False
                self._k_flash_arr = self._k_l_arr.copy()
                self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)
                logger.debug(f"Двухфазная (dew-наклон): S_v/S_l={ratio:.2f}")
            else:
                # Неоднозначная область → используем произведение (flash)
                self.stable = False
                self._k_flash_arr = self._k_v_arr * self._k_l_arr
                self.k_values_for_flash = self._array_to_dict(self._k_flash_arr)
                logger.debug(f"Неоднозначная область: используем flash K, S_v/S_l={ratio:.2f}")
            return
        
        # Fallback (не должно достигаться)
        logger.warning(f"Необычная ситуация: S_v={self.S_v:.3f}, S_l={self.S_l:.3f}")
        self.stable = False
        self.k_values_for_flash = self._array_to_dict(self._k_v_arr)

    # =====================================================================
    # ЦИКЛЫ ПО VAPOUR / LIQUID ТЕСТАМ
    # =====================================================================
    
    def _loop_vapour(self, k_values: np.ndarray, max_iter: int = 1000):
        self._k_v_arr = k_values.copy()
        
        for i in range(max_iter):
            Yi_v = self._calc_Yi_v(self._k_v_arr)
            self.S_v = self._calc_S(Yi_v)
            self._yi_v_arr = self._normalize_mole_fractions(Yi_v, self.S_v)
            self.yi_v = self._array_to_dict(self._yi_v_arr)
            
            new_comp = self._composition.new_composition(self.yi_v)
            eos = BrusilovskiyEOS(composition=new_comp, p=self.p, t=self.t)
            eos.calc_eos()
            
            ri_v = self._calc_ri_vapour(eos.fugacities, self._mixture_fugacities_arr, self.S_v)
            
            if self._check_convergence(ri_v):
                self.convergence_v = True
                self.k_values_vapour = self._array_to_dict(self._k_v_arr)
                return
            
            self._k_v_arr = self._update_k_values(self._k_v_arr, ri_v)
            
            if self._check_trivial_solution(self._k_v_arr):
                self.convergence_trivial_solution_v = True
                self.k_values_vapour = self._array_to_dict(self._k_v_arr)
                return
        
        logger.warning(f"Vapour-тест: достигнут лимит итераций ({max_iter})")
        self.k_values_vapour = self._array_to_dict(self._k_v_arr)
    
    def _loop_liquid(self, k_values: np.ndarray, max_iter: int = 1000):
        self._k_l_arr = k_values.copy()
        
        for i in range(max_iter):
            Xi_l = self._calc_Xi_l(self._k_l_arr)
            self.S_l = self._calc_S(Xi_l)
            self._xi_l_arr = self._normalize_mole_fractions(Xi_l, self.S_l)
            self.xi_l = self._array_to_dict(self._xi_l_arr)
            
            new_comp = self._composition.new_composition(self.xi_l)
            eos = BrusilovskiyEOS(composition=new_comp, p=self.p, t=self.t)
            eos.calc_eos()
            
            ri_l = self._calc_ri_liquid(eos.fugacities, self._mixture_fugacities_arr, self.S_l)
            
            if self._check_convergence(ri_l):
                self.convergence_l = True
                self.k_values_liquid = self._array_to_dict(self._k_l_arr)
                return
            
            self._k_l_arr = self._update_k_values(self._k_l_arr, ri_l)
            
            if self._check_trivial_solution(self._k_l_arr):
                self.convergence_trivial_solution_l = True
                self.k_values_liquid = self._array_to_dict(self._k_l_arr)
                return
        
        logger.warning(f"Liquid-тест: достигнут лимит итераций ({max_iter})")
        self.k_values_liquid = self._array_to_dict(self._k_l_arr)

    # =====================================================================
    # ОСНОВНОЙ РАСЧЕТ
    # =====================================================================
    
    def calculate_phase_stability(self):
        """Запускает полный тест стабильности."""
        # Фугитивности исходной смеси
        eos = BrusilovskiyEOS(composition=self._composition, p=self.p, t=self.t)
        eos.calc_eos()
        self._mixture_fugacities_arr = eos.fugacities.copy()
        
        # Начальные K по Вильсону
        k_init = self._calc_k_values_wilson_array()
        
        # Запускаем оба теста
        self._loop_vapour(k_init.copy())
        self._loop_liquid(k_init.copy())
        
        # Интерпретируем результаты
        self._interpretate_stability_analysis()

    # =====================================================================
    # СВОЙСТВО ДЛЯ РАСЧЕТА ТОЧКИ НАСЫЩЕНИЯ (УЛУЧШЕННОЕ)
    # =====================================================================
    
    @property
    def k_vals_for_sat_point_calculation(self) -> dict | None:
        """
        Возвращает K-факторы для расчета точки насыщения.
        
        Логика выбора:
        - Bubble point: S_v > 1, S_l <= 1 → k_values_vapour
        - Dew point: S_l > 1, S_v <= 1 → k_values_liquid
        - Неоднозначная область: по соотношению S_v/S_l
        """
        if (self.stable is None) or self.stable:
            return None
        
        # Чистый bubble point
        if (self.S_v > 1.0) and (self.S_l <= 1.0):
            logger.debug(f"✓ Bubble point K: S_v={self.S_v:.3f}, S_l={self.S_l:.3f}")
            return self.k_values_vapour
        
        # Чистый dew point
        if (self.S_l > 1.0) and (self.S_v <= 1.0):
            logger.debug(f"✓ Dew point K: S_l={self.S_l:.3f}, S_v={self.S_v:.3f}")
            return self.k_values_liquid
        
        # Оба > 1: выбираем по соотношению
        if (self.S_v > 1.0) and (self.S_l > 1.0):
            ratio = self.S_v / self.S_l if self.S_l > 1e-10 else float('inf')
            
            if ratio > 1.5:
                logger.debug(f"✓ Bubble-наклон: S_v/S_l={ratio:.2f} → vapour K")
                return self.k_values_vapour
            elif ratio < 0.67:
                logger.debug(f"✓ Dew-наклон: S_v/S_l={ratio:.2f} → liquid K")
                return self.k_values_liquid
            else:
                logger.debug(f"⚠ Неоднозначно: S_v/S_l={ratio:.2f} → flash K")
                return self.k_values_for_flash
        
        # Fallback
        logger.warning(f"⚠ Fallback: S_v={self.S_v:.3f}, S_l={self.S_l:.3f}")
        return self.k_values_for_flash