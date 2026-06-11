import numpy as np
import logging
from typing import Optional, Dict, Any, Tuple
from scipy.optimize import minimize, brentq
from _src.Composition.CompositionV2 import Composition
from _src.EOS.BrusilovskiyEOSV2 import BrusilovskiyEOS

logger = logging.getLogger(__name__)


def setup_logger(level: int = logging.INFO):
    """Настраивает логгер для вывода сообщений в консоль."""
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(level)
        formatter = logging.Formatter(
            '%(asctime)s | %(name)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(level)


class CriticalPointCalculator:
    """
    Расчет критической точки через термодинамические условия:
    (∂P/∂V)_T = 0 и (∂²P/∂V²)_T = 0
    
    Это классическое условие критической точки, не требующее расчета bubble/dew.
    """
    
    def __init__(self, composition: Composition, verbose: bool = True):
        if verbose:
            setup_logger(logging.INFO)
        
        self.composition = composition
        self.verbose = verbose
        
        # Результаты
        self.T_crit = None
        self.P_crit = None
        self.V_crit = None
        
    def _calc_pressure_from_eos(self, T: float, V: float) -> float:
        """
        Рассчитывает давление из EOS при заданных T и V.
        """
        try:
            eos = BrusilovskiyEOS(self.composition, 1.0, T)  # P не важен, пересчитаем
            
            # Безразмерные параметры
            A_m = eos._calc_A_mixture()
            B_m = eos._calc_B_mixture()
            C_m = eos._calc_C_mixture()
            D_m = eos._calc_D_mixture()
            
            # Z = PV/RT => P = ZRT/V
            # Кубическое уравнение: Z³ + E0*Z² + E1*Z + E2 = 0
            # где E0 = C + D - B - 1, E1 = A - BC + CD - BD - D - C, E2 = -(BCD + CD + AB)
            E_0 = C_m + D_m - B_m - 1.0
            E_1 = A_m - B_m * C_m + C_m * D_m - B_m * D_m - D_m - C_m
            E_2 = -(B_m * C_m * D_m + C_m * D_m + A_m * B_m)
            
            # Решаем для Z при заданном V
            # Z = PV/RT, но мы ищем P, поэтому используем итерационный подход
            # Или проще: используем явное выражение P(V,T) из EOS
            
            # Для уравнения Брусиловского:
            # P = RT/(V-b) - a/(V² + (b+c)V + bc)
            # где a, b, c, d - параметры смеси
            
            R = 83.14  # J/(mol·K) = cm³·bar/(mol·K)
            
            # Параметры смеси (упрощенно)
            a_mix = A_m * (R * T)**2 / 1.0  # нормализация
            b_mix = B_m * R * T / 1.0
            c_mix = C_m * R * T / 1.0
            d_mix = D_m * R * T / 1.0
            
            P = R * T / (V - b_mix) - a_mix / (V**2 + (b_mix + c_mix) * V + b_mix * c_mix)
            
            return P
            
        except Exception as e:
            logger.debug(f"EOS calculation failed at T={T}, V={V}: {e}")
            return np.inf
    
    def _calc_dP_dV(self, T: float, V: float, dV: float = 1e-6) -> float:
        """Первая производная dP/dV."""
        P1 = self._calc_pressure_from_eos(T, V - dV)
        P2 = self._calc_pressure_from_eos(T, V + dV)
        return (P2 - P1) / (2 * dV)
    
    def _calc_d2P_dV2(self, T: float, V: float, dV: float = 1e-6) -> float:
        """Вторая производная d²P/dV²."""
        dP1 = self._calc_dP_dV(T, V - dV, dV)
        dP2 = self._calc_dP_dV(T, V + dV, dV)
        return (dP2 - dP1) / (2 * dV)
    
    def _critical_conditions_objective(self, x: np.ndarray) -> float:
        """
        Целевая функция: сумма квадратов условий критической точки.
        x = [T, V]
        """
        T, V = x
        
        if T < 200 or T > 1000 or V < 0.01 or V > 10:
            return 1e10
        
        dP_dV = self._calc_dP_dV(T, V)
        d2P_dV2 = self._calc_d2P_dV2(T, V)
        
        # Нормализация для улучшения сходимости
        return (dP_dV / 100.0)**2 + (d2P_dV2 / 1000.0)**2
    
    def _scan_for_initial_guess(self, T_min: float, T_max: float, n_points: int = 30) -> Optional[Tuple[float, float]]:
        """
        Сканирование для нахождения начального приближения.
        Ищем температуру, где dP/dV ≈ 0 при некотором V.
        """
        logger.info(f"Scanning for initial guess in T=[{T_min-273.15:.1f}, {T_max-273.15:.1f}] °C")
        
        best_T = None
        best_V = None
        best_metric = float('inf')
        
        T_grid = np.linspace(T_min, T_max, n_points)
        
        for T in T_grid:
            # Для каждой температуры ищем V, где dP/dV ≈ 0
            V_grid = np.linspace(0.1, 5.0, 50)
            
            for V in V_grid:
                dP_dV = self._calc_dP_dV(T, V)
                d2P_dV2 = self._calc_d2P_dV2(T, V)
                
                metric = abs(dP_dV) + abs(d2P_dV2) * 0.1
                
                if metric < best_metric:
                    best_metric = metric
                    best_T = T
                    best_V = V
        
        if best_T is not None:
            logger.info(f"Best initial guess: T={best_T-273.15:.2f} °C, V={best_V:.4f}, metric={best_metric:.6f}")
            return best_T, best_V
        else:
            return None
    
    def calculate(self) -> Optional[Dict[str, Any]]:
        """
        Основной метод расчета критической точки.
        """
        logger.info("Starting critical point calculation via thermodynamic conditions...")
        
        # Определяем диапазон поиска
        comp_data = self.composition.composition_data
        z = np.array([self.composition.composition[c] for c in self.composition.composition.keys()])
        Tc_components = np.array([comp_data['critical_temperature'][c] for c in self.composition.composition.keys()])
        Tc_pseudo = np.sum(z * Tc_components)
        
        T_min = Tc_pseudo * 0.7
        T_max = Tc_pseudo * 1.5
        
        logger.info(f"Search range: T=[{T_min-273.15:.1f}, {T_max-273.15:.1f}] °C")
        
        # Шаг 1: Поиск начального приближения
        initial_guess = self._scan_for_initial_guess(T_min, T_max, n_points=40)
        
        if initial_guess is None:
            logger.error("Failed to find initial guess")
            return None
        
        T0, V0 = initial_guess
        
        # Шаг 2: Оптимизация
        logger.info(f"Starting optimization from T={T0-273.15:.2f} °C, V={V0:.4f}")
        
        try:
            result = minimize(
                self._critical_conditions_objective,
                x0=[T0, V0],
                method='Nelder-Mead',
                tol=1e-8,
                options={'maxiter': 1000, 'xatol': 1e-6, 'fatol': 1e-8}
            )
            
            if result.success or result.fun < 1e-6:
                T_crit, V_crit = result.x
                
                # Шаг 3: Расчет давления в критической точке
                P_crit = self._calc_pressure_from_eos(T_crit, V_crit)
                
                # Проверка условий
                dP_dV = self._calc_dP_dV(T_crit, V_crit)
                d2P_dV2 = self._calc_d2P_dV2(T_crit, V_crit)
                
                logger.info(f"✓ Critical point found:")
                logger.info(f"  T = {T_crit:.4f} K ({T_crit - 273.15:.2f} °C)")
                logger.info(f"  P = {P_crit:.4f} bar")
                logger.info(f"  V = {V_crit:.4f}")
                logger.info(f"  dP/dV = {dP_dV:.2e}")
                logger.info(f"  d²P/dV² = {d2P_dV2:.2e}")
                
                self.T_crit = T_crit
                self.P_crit = P_crit
                self.V_crit = V_crit
                
                return {
                    'T': T_crit,
                    'P': P_crit,
                    'V': V_crit,
                    'T_C': T_crit - 273.15,
                    'dP_dV': dP_dV,
                    'd2P_dV2': d2P_dV2,
                }
            else:
                logger.warning(f"Optimization did not converge well. Metric: {result.fun:.6f}")
                return None
                
        except Exception as e:
            logger.error(f"Optimization failed: {e}")
            return None
