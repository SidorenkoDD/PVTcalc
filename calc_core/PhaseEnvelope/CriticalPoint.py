import numpy as np
import logging
from typing import Optional, Dict, Any
from scipy.optimize import minimize_scalar
from calc_core.Composition.Composition import Composition
from calc_core.PhaseEnvelope.BubblePointPressure import BubblePointCalculator
from calc_core.PhaseEnvelope.DewPressure import DewPointCalculator

logger = logging.getLogger(__name__)


class CriticalPointCalculator:
    """
    Расчет критической точки через минимизацию зазора между bubble и dew кривыми.
    """

    def __init__(self, composition: Composition, verbose: bool = True):
        # Настройка вывода логов (уровень/handler) — забота вызывающего кода,
        # см. calc_core/Utils/Logging.py::configure_logging(). verbose оставлен
        # только для обратной совместимости сигнатуры, сам по себе больше
        # ничего не настраивает.
        self.composition = composition
        self.verbose = verbose
        self.T_crit = None
        self.P_crit = None
        
    def _calc_gap(self, T: float) -> float:
        """
        Вычисляет |P_dew(T) - P_bubble(T)|.
        Возвращает 1e6, если расчет некорректен.
        """
        try:
            self.composition.T = T
            
            # Bubble point
            b_calc = BubblePointCalculator(
                composition=self.composition,
                T=T,
                max_iter=150,
                tol=1e-8
            )
            P_b = b_calc.calculate()
            if not b_calc.converged or P_b is None or P_b <= 0:
                return 1e6
            
            # Dew point (без dew_point_type, так как его нет в текущей версии)
            d_calc = DewPointCalculator(
                composition=self.composition,
                T=T,
                max_iter=150,
                tol=1e-8
            )
            P_d = d_calc.calculate()
            
            if P_d is None or P_d <= 0:
                return 1e6
            
            # Проверка: P_dew должен быть >= P_bubble
            if P_d < P_b * 0.9:
                return 1e6
            
            return abs(P_d - P_b)
            
        except Exception:
            return 1e6
    
    def calculate(self) -> Optional[Dict[str, Any]]:
        """
        Основной метод расчета критической точки.
        """
        logger.info("Starting critical point calculation via bubble-dew gap minimization...")
        
        # Определяем диапазон поиска
        comp_data = self.composition.composition_data
        z = np.array([self.composition.composition[c] for c in self.composition.composition.keys()])
        Tc_comps = np.array([comp_data['critical_temperature'][c] for c in self.composition.composition.keys()])
        Tc_pseudo = np.sum(z * Tc_comps)
        
        T_min = Tc_pseudo * 0.85
        T_max = Tc_pseudo * 1.25
        
        logger.info(f"Search range: T=[{T_min-273.15:.1f}, {T_max-273.15:.1f}] °C")
        
        # Шаг 1: Грубое сканирование с шагом 2K
        T_grid = np.arange(T_min, T_max, 2.0)
        best_T = None
        best_gap = 1e6
        
        logger.info(f"Scanning {len(T_grid)} points...")
        
        for T in T_grid:
            gap = self._calc_gap(T)
            if gap < best_gap:
                best_gap = gap
                best_T = T
        
        if best_T is None or best_gap > 1e5:
            logger.error("Failed to find critical region")
            return None
        
        logger.info(f"Best guess: T={best_T-273.15:.2f} °C, gap={best_gap:.4f} bar")
        
        # Шаг 2: Уточнение методом Brent на интервале [T_best - 5, T_best + 5]
        logger.info("Refining critical point...")
        try:
            res = minimize_scalar(
                self._calc_gap,
                bounds=(best_T - 5.0, best_T + 5.0),
                method='bounded',
                options={'xatol': 1e-4, 'maxiter': 50}
            )
            
            if res.success and res.fun < 1e5:
                T_refined = res.x
            else:
                T_refined = best_T
        except Exception:
            T_refined = best_T
        
        # Шаг 3: Финальный расчет давлений
        self.composition.T = T_refined
        
        b_calc = BubblePointCalculator(
            composition=self.composition,
            T=T_refined,
            max_iter=150,
            tol=1e-8
        )
        P_b = b_calc.calculate()
        
        d_calc = DewPointCalculator(
            composition=self.composition,
            T=T_refined,
            max_iter=150,
            tol=1e-8
        )
        P_d = d_calc.calculate()
        
        if P_b is None or P_d is None or P_b <= 0 or P_d <= 0:
            logger.error("Final calculation failed")
            return None
        
        P_crit = (P_b + P_d) / 2.0
        gap_final = abs(P_b - P_d)
        
        logger.info(f"✓ Critical point found:")
        logger.info(f"  T = {T_refined:.4f} K ({T_refined - 273.15:.2f} °C)")
        logger.info(f"  P = {P_crit:.4f} bar")
        logger.info(f"  P_bubble = {P_b:.4f} bar")
        logger.info(f"  P_dew = {P_d:.4f} bar")
        logger.info(f"  Gap = {gap_final:.6f} bar")
        
        self.T_crit = T_refined
        self.P_crit = P_crit
        
        return {
            'T': T_refined,
            'P': P_crit,
            'T_C': T_refined - 273.15,
            'P_bubble': P_b,
            'P_dew': P_d,
            'gap': gap_final
        }


if __name__ == "__main__":
    from calc_core.Composition.Composition import Composition
    
    composition_dict = {
        'N2': 0.01, 'C1': 0.70, 'C2': 0.08, 'C3': 0.05,
        'iC4': 0.01, 'nC4': 0.02, 'iC5': 0.01, 'nC5': 0.01,
        'C6': 0.01, 'C7': 0.05, 'C10': 0.03, 'C15': 0.01,
        'C20': 0.005, 'C30': 0.005
    }
    
    composition = Composition(composition_dict, T_res=350.0)
    composition.evaluate_composition_data(c7_plus_correlations={})
    
    crit_calc = CriticalPointCalculator(composition, verbose=True)
    crit_point = crit_calc.calculate()
    
    if crit_point:
        print(f"\n{'='*40}")
        print(f"CRITICAL POINT RESULTS")
        print(f"{'='*40}")
        print(f"T_crit = {crit_point['T_C']:.2f} °C")
        print(f"P_crit = {crit_point['P']:.2f} bar")
        print(f"Gap    = {crit_point['gap']:.4f} bar")
    else:
        print("\nCritical point calculation failed.")