import numpy as np
import logging
from _src.EOS.BrusilovskiyEOSV2 import BrusilovskiyEOS
from _src.Composition.CompositionV2 import Composition
from _src.PhaseStability.TwoPhaseStabilityTestV3 import TwoPhaseStabilityTest

logger = logging.getLogger(__name__)


class DewPointCalculator:
    """
    Класс для расчета давления начала конденсации (Dew Point Pressure).
    
    dew_point_type:
      'upper' - верхнее давление начала конденсации (существует при T > T_crit, высокое P)
      'lower' - нижнее давление начала конденсации (существует при T > T_crit, низкое P)
    """

    def __init__(
        self,
        composition: Composition,
        T: float,
        dew_point_type: str = 'upper',
        P_guess: float = None,
        max_iter: int = 200,
        tol: float = 1e-13,
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

        self.result_P = None
        self.result_F = None
        self.converged = False
        self.iterations_done = 0

    def calculate(self):
        """
        Выполняет расчет давления начала конденсации.
        """
        comp_data = self.composition.composition_data
        components = list(self.composition.composition.keys())
        z = np.array([self.composition.composition[c] for c in components])

        Pc = np.array([comp_data['critical_pressure'][c] for c in components])
        Tc = np.array([comp_data['critical_temperature'][c] for c in components])
        omega = np.array([comp_data['acentric_factor'][c] for c in components])

        exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))

        # Начальное приближение зависит от типа dew point
        if self.P_guess is not None:
            P = self.P_guess
        elif self.dew_point_type == 'upper':
            # Upper dew point: высокое начальное P
            P = np.sum(z * Pc * exp_term) * 1.5
        else:
            # Lower dew point: низкое начальное P
            K_wilson = (Pc / 1.0) * exp_term
            P = 1.0 / np.sum(z / K_wilson)

        K = (Pc / P) * exp_term

        # Находим стабильную область для начала итераций
        result = self._find_stable_region(P, K, z, components)
        if result is None:
            logger.info(f"{self.dew_point_type} dew point does not exist at T={self.T:.2f} K")
            self.result_P = None
            self.converged = False
            return None
        
        P, K = result

        # Основной цикл итераций
        F = None
        trivial_counter = 0
        
        for j in range(self.max_iter):
            # Состав жидкой фазы: x_i = z_i / K_i
            x = z / K
            x_sum = np.sum(x)
            x_normalized = x / x_sum
            x_dict = {c: float(x_normalized[i]) for i, c in enumerate(components)}
            liquid_comp = self.composition.new_composition(x_dict)
            liquid_comp.T = self.T

            # Расчет УРС для обеих фаз
            eos_V = BrusilovskiyEOS(self.composition, P, self.T)
            eos_L = BrusilovskiyEOS(liquid_comp, P, self.T)

            eos_V.calc_eos()
            eos_L.calc_eos()

            Z_V = np.max(eos_V.real_roots_eos)
            Z_L = np.min(eos_L.real_roots_eos)

            # Проверка на тривиальное решение
            if np.isclose(Z_L, Z_V, rtol=1e-4):
                trivial_counter += 1
                logger.warning(f"Iter {j}: Trivial solution detected (Z_L ≈ Z_V = {Z_L:.6f})")
                
                if trivial_counter >= 3:
                    logger.info(f"Iter {j}: Escaping trivial solution via stability test")
                    result = self._escape_trivial_solution(P, K, z, components)
                    
                    if result is None:
                        logger.warning("Cannot escape trivial solution. Returning None.")
                        self.result_P = None
                        self.converged = False
                        return None
                    
                    P_new, K_new = result
                    P = P_new
                    if K_new is not None:
                        K = K_new
                    trivial_counter = 0
                    continue

            ln_phi_V = eos_V.get_fugacity_coef_vector_by_root(Z_V)
            ln_phi_L = eos_L.get_fugacity_coef_vector_by_root(Z_L)

            # Новые K-факторы
            ln_K_new = ln_phi_L - ln_phi_V
            K_new = np.exp(ln_K_new)

            # Функция F
            F = np.sum(z / K_new) - 1.0

            if abs(F) < self.tol:
                logger.info(f"{self.dew_point_type} dew point converged in {j+1} iterations. P = {P:.4f}, F = {F:.2e}")
                self.result_P = P
                self.result_F = F
                self.converged = True
                self.iterations_done = j + 1
                return P

            # Производная dF/dP
            eos_V._z_factor = Z_V
            eos_V._dz_dp = None
            eos_V._dlogphi_dp = None
            dlnphiV_dP = np.array([eos_V.calc_d_log_phi_i_dp(c) for c in components])

            eos_L._z_factor = Z_L
            eos_L._dz_dp = None
            eos_L._dlogphi_dp = None
            dlnphiL_dP = np.array([eos_L.calc_d_log_phi_i_dp(c) for c in components])

            dF_dP = np.sum((z / K_new) * (dlnphiV_dP - dlnphiL_dP))

            if dF_dP < 0:
                logger.debug(f"Iter {j}: Sign correction for dF/dP. Original: {dF_dP:.2e} -> {abs(dF_dP):.2e}")
                dF_dP = abs(dF_dP)

            if abs(dF_dP) < 1e-15:
                logger.warning(f"dF/dP too small at iteration {j}. Stopping.")
                break

            # Обновление давления
            delta_P = -F / dF_dP
            max_step = 0.5 * P
            if abs(delta_P) > max_step:
                delta_P = max_step * np.sign(delta_P)

            P_new = P + delta_P
            if P_new <= 0:
                P_new = P * 0.5

            P = P_new
            K = K_new

            if j % 5 == 0:
                logger.debug(f"Iter {j}: P = {P:.4f}, F = {F:.6f}, dF/dP = {dF_dP:.6e}, Z_L = {Z_L:.6f}, Z_V = {Z_V:.6f}")

        logger.warning(f"Max iterations ({self.max_iter}) reached for {self.dew_point_type} dew point. P = {P:.4f}, F = {F:.2e}. Returning None")
        self.result_P = None
        self.result_F = F
        self.converged = False
        self.iterations_done = self.max_iter
        return None

    def _find_stable_region(self, P: float, K: np.ndarray, z: np.ndarray, components: list):
        """
        Находит стабильную область для начала итераций.
        """
        if self.dew_point_type == 'upper':
            return self._find_upper_dew_point_region(P, K, z, components)
        else:
            return self._find_lower_dew_point_region(P, K, z, components)

    def _find_upper_dew_point_region(self, P: float, K: np.ndarray, z: np.ndarray, components: list):
        """
        Поиск upper dew point: начинаем с высокого P и идём ВНИЗ до нестабильности.
        Upper dew point - это граница между стабильной жидкостью/сверхкритической фазой (высокое P)
        и двухфазной областью.
        """
        P_high = P  # Стабильная область (высокое P)
        P_low = P_high / 10.0
        
        # Проверяем стабильность при P_high
        try:
            st = TwoPhaseStabilityTest(self.composition, P_high, self.T)
            st.calculate_phase_stability()
        except Exception as e:
            logger.warning(f"Stability test failed at P={P_high}: {e}")
            return None
        
        if not st.stable:
            # Система нестабильна даже при высоком P - нужно идти еще выше
            logger.info(f"System unstable at P={P_high}. Searching higher...")
            for attempt in range(10):
                P_high *= 2.0
                if P_high > 5000.0:
                    logger.warning(f"Pressure exceeded 5000 bar. Upper dew point may not exist.")
                    return None
                try:
                    st = TwoPhaseStabilityTest(self.composition, P_high, self.T)
                    st.calculate_phase_stability()
                    if st.stable:
                        break
                except Exception as e:
                    return None
            
            if not st.stable:
                return None
        
        # Теперь P_high - стабильная область. Идём ВНИЗ пока не найдём нестабильность.
        unstable_found = False
        
        for attempt in range(20):
            if P_low < 0.01:
                logger.warning(f"Pressure below 0.01 bar. Upper dew point may not exist.")
                break
            
            try:
                st = TwoPhaseStabilityTest(self.composition, P_low, self.T)
                st.calculate_phase_stability()
                
                if not st.stable:
                    unstable_found = True
                    logger.info(f"Found unstable region at P={P_low:.4f}")
                    break
            except Exception as e:
                logger.warning(f"Stability test failed at P={P_low}: {e}")
                return None
            
            P_high = P_low
            P_low /= 2.0
        
        if not unstable_found:
            logger.info(f"System stable down to P={P_low}. Upper dew point does not exist at T={self.T:.2f} K")
            return None
        
        # Бисекция между нестабильным (P_low) и стабильным (P_high)
        K_stable = None
        
        for i in range(30):
            P_mid = (P_low + P_high) / 2.0
            
            try:
                st = TwoPhaseStabilityTest(self.composition, P_mid, self.T)
                st.calculate_phase_stability()
                
                if st.stable:
                    P_high = P_mid
                    if st.k_values_for_flash is not None:
                        K_stable = np.array([st.k_values_for_flash[c] for c in components])
                else:
                    P_low = P_mid
                
                if (P_high - P_low) / P_high < 0.05:
                    break
            except Exception as e:
                logger.warning(f"Stability test failed during bisection at P={P_mid}: {e}")
                break
        
        P_new = P_high  # Стабильное давление (выше upper dew point)
        
        if K_stable is not None:
            logger.info(f"Found stable region at P={P_new:.4f}. Using stability test K-values.")
            return P_new, K_stable
        else:
            logger.info(f"Found stable region at P={P_new:.4f}. Using Wilson K-values.")
            comp_data = self.composition.composition_data
            Pc = np.array([comp_data['critical_pressure'][c] for c in components])
            Tc = np.array([comp_data['critical_temperature'][c] for c in components])
            omega = np.array([comp_data['acentric_factor'][c] for c in components])
            exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
            K_new = (Pc / P_new) * exp_term
            return P_new, K_new

    def _find_lower_dew_point_region(self, P: float, K: np.ndarray, z: np.ndarray, components: list):
        """
        Поиск lower dew point: начинаем с низкого P и идём ВВЕРХ до нестабильности.
        Lower dew point - это граница между стабильным газом (низкое P) и двухфазной областью.
        """
        P_low = 0.01  # Начальное низкое давление (стабильный газ)
        P_high = P_low * 10.0
        
        # Проверяем стабильность при P_low
        try:
            st = TwoPhaseStabilityTest(self.composition, P_low, self.T)
            st.calculate_phase_stability()
        except Exception as e:
            logger.warning(f"Stability test failed at P={P_low}: {e}")
            return None
        
        if not st.stable:
            # Система нестабильна даже при очень низком P
            logger.info(f"System unstable at P={P_low}. Lower dew point does not exist at T={self.T:.2f} K")
            return None
        
        # Идём ВВЕРХ по давлению пока не найдём нестабильность
        unstable_found = False
        
        for attempt in range(20):
            if P_high > 1000.0:
                logger.warning(f"Pressure exceeded 1000 bar. Lower dew point may not exist.")
                break
            
            try:
                st = TwoPhaseStabilityTest(self.composition, P_high, self.T)
                st.calculate_phase_stability()
                
                if not st.stable:
                    unstable_found = True
                    logger.info(f"Found unstable region at P={P_high:.4f}")
                    break
            except Exception as e:
                logger.warning(f"Stability test failed at P={P_high}: {e}")
                return None
            
            P_low = P_high
            P_high *= 2.0
        
        if not unstable_found:
            logger.info(f"System stable up to P={P_high}. Lower dew point does not exist at T={self.T:.2f} K")
            return None
        
        # Бисекция между стабильным (P_low) и нестабильным (P_high)
        K_stable = None
        
        for i in range(30):
            P_mid = (P_low + P_high) / 2.0
            
            try:
                st = TwoPhaseStabilityTest(self.composition, P_mid, self.T)
                st.calculate_phase_stability()
                
                if st.stable:
                    P_low = P_mid
                    if st.k_values_for_flash is not None:
                        K_stable = np.array([st.k_values_for_flash[c] for c in components])
                else:
                    P_high = P_mid
                
                if (P_high - P_low) / P_high < 0.05:
                    break
            except Exception as e:
                logger.warning(f"Stability test failed during bisection at P={P_mid}: {e}")
                break
        
        P_new = P_low  # Стабильное давление (ниже lower dew point)
        
        if K_stable is not None:
            logger.info(f"Found stable region at P={P_new:.4f}. Using stability test K-values.")
            return P_new, K_stable
        else:
            logger.info(f"Found stable region at P={P_new:.4f}. Using Wilson K-values.")
            comp_data = self.composition.composition_data
            Pc = np.array([comp_data['critical_pressure'][c] for c in components])
            Tc = np.array([comp_data['critical_temperature'][c] for c in components])
            omega = np.array([comp_data['acentric_factor'][c] for c in components])
            exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
            K_new = (Pc / P_new) * exp_term
            return P_new, K_new

    def _escape_trivial_solution(self, P: float, K: np.ndarray, z: np.ndarray, components: list):
        """
        Пытается выйти из тривиального решения.
        """
        try:
            stability_test = TwoPhaseStabilityTest(self.composition, P, self.T)
            stability_test.calculate_phase_stability()
        except Exception as e:
            logger.warning(f"Stability test failed: {e}")
            return None
        
        if stability_test.stable:
            if self.dew_point_type == 'upper':
                P_new = P * 0.85
            else:
                P_new = P * 1.15
            
            logger.info(f"Stable at P={P:.4f}, adjusting to P={P_new:.4f}")
            
            comp_data = self.composition.composition_data
            Pc = np.array([comp_data['critical_pressure'][c] for c in components])
            Tc = np.array([comp_data['critical_temperature'][c] for c in components])
            omega = np.array([comp_data['acentric_factor'][c] for c in components])
            exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
            K_new = (Pc / P_new) * exp_term
            
            return P_new, K_new
        else:
            if stability_test.k_values_for_flash is not None:
                K_new = np.array([stability_test.k_values_for_flash[c] for c in components])
                return P, K_new
            else:
                P_new = P * 0.85
                comp_data = self.composition.composition_data
                Pc = np.array([comp_data['critical_pressure'][c] for c in components])
                Tc = np.array([comp_data['critical_temperature'][c] for c in components])
                omega = np.array([comp_data['acentric_factor'][c] for c in components])
                exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
                K_new = (Pc / P_new) * exp_term
                return P_new, K_new
