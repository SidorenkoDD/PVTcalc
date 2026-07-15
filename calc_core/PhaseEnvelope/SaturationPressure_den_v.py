import numpy as np
import logging
from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.Composition.Composition import Composition
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest

logger = logging.getLogger(__name__)


def norm_array(x: np.ndarray) -> np.ndarray:
    """Нормализация массива"""
    s = float(np.sum(x))
    return x / s


def Y_from_K_dew_array(z: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Расчет состава жидкой фазы для dew point: x_i = z_i / K_i"""
    return z / K


def check_convergence_array(z: np.ndarray, Y: np.ndarray, R: np.ndarray) -> bool:
    """Проверка сходимости"""
    flag1 = abs(1.0 - float(np.sum(Y))) < 1e-10
    ratio = Y / z
    
    # Защита от деления на ноль
    denom = np.log(ratio)
    numer = np.log(R)
    
    # Фильтруем малые значения
    valid = np.abs(denom) > 1e-30
    if not np.any(valid):
        flag2 = True
    else:
        metric = float(np.sum(numer[valid] / denom[valid]))
        flag2 = (metric ** 2) < 1e-6
    
    return flag1 and flag2


def check_trivial_array(z: np.ndarray, Y: np.ndarray) -> bool:
    """Проверка на тривиальное решение"""
    metric = float(np.sum(np.log(Y / z)))
    return (metric ** 2) < 1e-4


class DewPointCalculator:
    """
    Класс для расчета давления начала конденсации (Dew Point Pressure)
    по алгоритму из SaturationPressure_den_v.py
    
    dew_point_type:
      'upper' - верхнее давление начала конденсации
      'lower' - нижнее давление начала конденсации
    """

    def __init__(
        self,
        composition: Composition,
        T: float,
        dew_point_type: str = 'upper',
        P_guess: float = None,
        max_iter: int = 200,
        tol: float = 1e-10,
        damp: float = 0.6,
    ):
        if dew_point_type not in ['upper', 'lower']:
            raise ValueError(f"dew_point_type must be 'upper' or 'lower', got '{dew_point_type}'")
        
        self.composition = composition
        self.T = T
        self.composition.T = self.T
        self.dew_point_type = dew_point_type
        self.P_guess = P_guess
        self.max_iter = max_iter
        self.tol = tol
        self.damp = damp  # Коэффициент демпфирования

        self.result_P = None
        self.converged = False
        self.iterations_done = 0

    def calculate(self) -> float:
        """
        Расчет давления начала конденсации.
        Использует алгоритм с bracketing (p1, p2) и итерациями по составу Y.
        """
        logger.info(f"Starting dew point calculation at T={self.T:.2f} K, type={self.dew_point_type}")
        
        # Получаем данные компонент
        comp_data = self.composition.composition_data
        components = list(self.composition.composition.keys())
        z = np.array([self.composition.composition[c] for c in components])
        
        # Критические параметры для начального приближения
        Pc = np.array([comp_data['critical_pressure'][c] for c in components])
        Tc = np.array([comp_data['critical_temperature'][c] for c in components])
        omega = np.array([comp_data['acentric_factor'][c] for c in components])
        
        # Шаг 1: Находим границы p1 (стабильное) и p2 (нестабильное)
        p1, p2, K_init = self._find_pressure_bounds(z, components, Pc, Tc, omega)
        
        if p1 is None or p2 is None:
            logger.warning(f"Cannot find pressure bounds for {self.dew_point_type} dew point")
            self.result_P = None
            self.converged = False
            return None
        
        logger.info(f"Pressure bounds: p1 (stable)={p1:.4f} bar, p2 (unstable)={p2:.4f} bar")
        
        # Шаг 2: Инициализация состава Y из K-факторов
        K_arr = K_init
        Y = Y_from_K_dew_array(z, K_arr)
        Y = norm_array(Y)
        
        logger.debug(f"Initial Y composition: {dict(zip(components, Y))}")
        
        # Начальное давление (между p1 и p2)
        p = p2  # Начинаем с нестабильной области
        
        # Шаг 3: Основной итерационный процесс
        for i in range(self.max_iter):
            logger.debug(f"\n=== Iteration {i+1}/{self.max_iter} ===")
            logger.debug(f"Current pressure: p={p:.6f} bar")
            
            # Нормализуем Y
            S = float(np.sum(Y))
            y = norm_array(Y)
            y_dict = {c: float(y[j]) for j, c in enumerate(components)}
            
            logger.debug(f"Sum Y: S={S:.6f}")
            
            # EOS для исходной смеси z (паровая фаза)
            eos_z = BrusilovskiyEOS(self.composition, p, self.T)
            eos_z.calc_eos()
            Z_z = np.max(eos_z.real_roots_eos)
            ln_phi_z = eos_z.get_fugacity_coef_vector_by_root(Z_z)
            
            # EOS для зарождающейся жидкой фазы y
            y_composition = self.composition.new_composition(y_dict)
            y_composition.T = self.T
            eos_y = BrusilovskiyEOS(y_composition, p, self.T)
            eos_y.calc_eos()
            Z_y = np.min(eos_y.real_roots_eos)
            ln_phi_y = eos_y.get_fugacity_coef_vector_by_root(Z_y)
            
            logger.debug(f"Z_z={Z_z:.6f}, Z_y={Z_y:.6f}")
            logger.debug(f"ln_phi_z: {dict(zip(components, ln_phi_z))}")
            logger.debug(f"ln_phi_y: {dict(zip(components, ln_phi_y))}")
            
            # Проверка на тривиальное решение
            if np.isclose(Z_z, Z_y, rtol=1e-4):
                logger.warning(f"Trivial solution detected at iteration {i+1}")
                p, Y, K_arr = self._escape_trivial(p, Y, z, K_arr, components, Pc, Tc, omega)
                if p is None:
                    return None
                continue
            
            # Расчет R_i = (z_i * φ_i(z)) / (y_i * φ_i(y)) / S
            # Работаем в логарифмах: ln(R_i) = ln(z_i) + ln_phi_z - ln(y_i) - ln_phi_y - ln(S)
            R = np.exp(
                np.log(z) + ln_phi_z - np.log(y) - ln_phi_y
            ) / S
            
            logger.debug(f"R factors: {dict(zip(components, R))}")
            
            # Обновление Y: Y_new = Y * R
            Y_new = Y * R
            S_new = float(np.sum(Y_new))
            
            # Проверка сходимости
            convergence = check_convergence_array(z, Y_new, R)
            trivial = check_trivial_array(z, Y_new)
            
            logger.info(f"Iteration {i+1}: p={p:.6f} bar, S={S_new:.6f}, convergence={convergence}")
            
            if convergence:
                logger.info(f"Dew point converged in {i+1} iterations. P={p:.4f} bar")
                self.result_P = p
                self.converged = True
                self.iterations_done = i + 1
                return p
            
            if trivial:
                logger.warning("Trivial solution detected. Adjusting pressure.")
                p, Y, K_arr = self._escape_trivial(p, Y, z, K_arr, components, Pc, Tc, omega)
                if p is None:
                    return None
                continue
            
            # Расчет производной dQ/dp
            # Q = 1.0 - S
            # dQdp = -dS/dp
            dlnphi_z_dp = eos_z._calc_dlogphi_dp_vector()
            dlnphi_y_dp = eos_y._calc_dlogphi_dp_vector()
            
            # dQdp = sum((Y * R) * (dlnphi_y_dp - dlnphi_z_dp))
            dQdp = float(np.sum((Y * R) * (dlnphi_y_dp - dlnphi_z_dp)))
            
            logger.debug(f"dQdp={dQdp:.6e}")
            
            Q = 1.0 - S
            
            # Метод Ньютона для обновления давления
            if abs(dQdp) > 1e-20 and np.isfinite(dQdp):
                delta_p = -self.damp * Q / dQdp
                
                # Ограничиваем шаг
                max_step = 0.3 * p
                if abs(delta_p) > max_step:
                    delta_p = max_step * np.sign(delta_p)
                
                p_new = p + delta_p
                
                # Проверяем границы
                if p_new < p2:
                    p_new = p2 * 1.01
                if p_new > p1:
                    p_new = p1 * 0.99
                
                if p_new <= 0:
                    p_new = p * 0.5
                
                logger.debug(f"Newton step: Q={Q:.6f}, dQdp={dQdp:.6e}, delta_p={delta_p:.6f}")
                logger.debug(f"Pressure update: {p:.6f} -> {p_new:.6f}")
                
                p = p_new
            else:
                logger.warning(f"dQdp too small or infinite. Using bisection.")
                # Fallback к бисекции
                p = (p1 + p2) / 2.0
            
            # Обновляем Y
            Y = norm_array(Y_new)
            
            # Обновляем границы если нужно
            if S > 1.0:
                # Давление слишком низкое (нестабильная область)
                p2 = max(p2, p)
            else:
                # Давление слишком высокое (стабильная область)
                p1 = min(p1, p)
        
        logger.warning(f"Max iterations ({self.max_iter}) reached. P={p:.4f} bar")
        self.result_P = None
        self.converged = False
        self.iterations_done = self.max_iter
        return None

    def _find_pressure_bounds(self, z, components, Pc, Tc, omega):
        """
        Находит границы p1 (стабильное) и p2 (нестабильное) для dew point.
        """
        exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
        
        if self.dew_point_type == 'upper':
            # Upper dew point: высокое давление
            # p1 - стабильное (высокое P), p2 - нестабильное (ниже)
            
            # Начальное приближение
            if self.P_guess is not None:
                P_init = self.P_guess
            else:
                P_init = np.sum(z * Pc * exp_term) * 1.5
            
            # Проверяем стабильность при P_init
            try:
                st = TwoPhaseStabilityTest(self.composition, P_init, self.T)
                st.calculate_phase_stability()
                
                if st.stable:
                    # P_init стабильно - это p1
                    p1 = P_init
                    # Ищем p2 (нестабильное)
                    p2 = p1 / 5.0
                    for _ in range(10):
                        st2 = TwoPhaseStabilityTest(self.composition, p2, self.T)
                        st2.calculate_phase_stability()
                        if not st2.stable:
                            break
                        p1 = p2
                        p2 /= 2.0
                    else:
                        logger.warning("Cannot find unstable region for upper dew point")
                        return None, None, None
                else:
                    # P_init нестабильно - это p2
                    p2 = P_init
                    # Ищем p1 (стабильное)
                    p1 = p2 * 5.0
                    for _ in range(10):
                        st1 = TwoPhaseStabilityTest(self.composition, p1, self.T)
                        st1.calculate_phase_stability()
                        if st1.stable:
                            break
                        p2 = p1
                        p1 *= 2.0
                    else:
                        logger.warning("Cannot find stable region for upper dew point")
                        return None, None, None
                
                # Получаем K-факторы из теста стабильности
                if st.k_values_for_flash is not None:
                    K_init = np.array([st.k_values_for_flash[c] for c in components])
                else:
                    K_init = (Pc / p2) * exp_term
                
                logger.info(f"Upper dew point bounds: p1={p1:.4f}, p2={p2:.4f}")
                return p1, p2, K_init
                
            except Exception as e:
                logger.error(f"Error finding pressure bounds: {e}")
                return None, None, None
                
        else:
            # Lower dew point: низкое давление
            # p1 - стабильное (низкое P), p2 - нестабильное (выше)
            
            if self.P_guess is not None:
                P_init = self.P_guess
            else:
                K_wilson = (Pc / 1.0) * exp_term
                P_init = 1.0 / np.sum(z / K_wilson)
            
            try:
                st = TwoPhaseStabilityTest(self.composition, P_init, self.T)
                st.calculate_phase_stability()
                
                if st.stable:
                    # P_init стабильно - это p1
                    p1 = P_init
                    # Ищем p2 (нестабильное)
                    p2 = p1 * 2.0
                    for _ in range(10):
                        st2 = TwoPhaseStabilityTest(self.composition, p2, self.T)
                        st2.calculate_phase_stability()
                        if not st2.stable:
                            break
                        p1 = p2
                        p2 *= 2.0
                    else:
                        logger.warning("Cannot find unstable region for lower dew point")
                        return None, None, None
                else:
                    # P_init нестабильно - это p2
                    p2 = P_init
                    # Ищем p1 (стабильное)
                    p1 = p2 / 5.0
                    for _ in range(10):
                        st1 = TwoPhaseStabilityTest(self.composition, p1, self.T)
                        st1.calculate_phase_stability()
                        if st1.stable:
                            break
                        p2 = p1
                        p1 /= 2.0
                    else:
                        logger.warning("Cannot find stable region for lower dew point")
                        return None, None, None
                
                if st.k_values_for_flash is not None:
                    K_init = np.array([st.k_values_for_flash[c] for c in components])
                else:
                    K_init = (Pc / p1) * exp_term
                
                logger.info(f"Lower dew point bounds: p1={p1:.4f}, p2={p2:.4f}")
                return p1, p2, K_init
                
            except Exception as e:
                logger.error(f"Error finding pressure bounds: {e}")
                return None, None, None

    def _escape_trivial(self, p, Y, z, K_arr, components, Pc, Tc, omega):
        """Выход из тривиального решения"""
        try:
            st = TwoPhaseStabilityTest(self.composition, p, self.T)
            st.calculate_phase_stability()
            
            if st.k_values_for_flash is not None:
                K_new = np.array([st.k_values_for_flash[c] for c in components])
                Y_new = Y_from_K_dew_array(z, K_new)
                Y_new = norm_array(Y_new)
                
                if self.dew_point_type == 'upper':
                    p_new = p * 0.9
                else:
                    p_new = p * 1.1
                
                logger.info(f"Escaping trivial solution using stability K. P: {p:.4f} -> {p_new:.4f}")
                return p_new, Y_new, K_new
        except Exception as e:
            logger.warning(f"Stability test failed in escape: {e}")
        
        # Fallback
        exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
        if self.dew_point_type == 'upper':
            p_new = p * 0.85
        else:
            p_new = p * 1.15
        
        K_new = (Pc / p_new) * exp_term
        Y_new = Y_from_K_dew_array(z, K_new)
        Y_new = norm_array(Y_new)
        
        logger.info(f"Escaping trivial solution (fallback). P: {p:.4f} -> {p_new:.4f}")
        return p_new, Y_new, K_new


# Обратная совместимость
def calculate_dew_point_pressure(
    composition: Composition,
    T: float,
    dew_point_type: str = 'upper',
    P_guess: float = None,
    max_iter: int = 200,
    tol: float = 1e-10,
) -> float:
    calculator = DewPointCalculator(
        composition=composition,
        T=T,
        dew_point_type=dew_point_type,
        P_guess=P_guess,
        max_iter=max_iter,
        tol=tol,
    )
    return calculator.calculate()