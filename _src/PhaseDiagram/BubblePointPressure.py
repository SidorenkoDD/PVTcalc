import numpy as np
import logging
from _src.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from _src.Composition.Composition import Composition
from _src.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest

logger = logging.getLogger(__name__)


class BubblePointCalculator:
    """
    Класс для расчета давления точки кипения (Bubble Point Pressure)
    по алгоритму из Table 6.1 с интеграцией теста стабильности.
    """

    
    def __init__(
        self,
        composition: Composition,
        T: float,
        P_guess: float = None,
        max_iter: int = 200,
        tol: float = 1e-8,
    ):
        self.composition = composition
        self.T = T
        self.composition.T = self.T
        self.P_guess = P_guess
        self.max_iter = max_iter
        self.tol = tol

        self.result_P = None
        self.result_F = None
        self.converged = False
        self.iterations_done = 0

    def calculate(self) -> float:
        """
        Выполняет расчет давления точки кипения.
        """
        comp_data = self.composition.composition_data
        components = list(self.composition.composition.keys())
        z = np.array([self.composition.composition[c] for c in components])

        Pc = np.array([comp_data['critical_pressure'][c] for c in components])
        Tc = np.array([comp_data['critical_temperature'][c] for c in components])
        omega = np.array([comp_data['acentric_factor'][c] for c in components])

        # 1. Начальное приближение (Wilson)
        exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))

        if self.P_guess is None:
            P = np.sum(z * Pc * exp_term)
        else:
            P = self.P_guess

        K = (Pc / P) * exp_term

        # Проверяем стабильность при начальном P и находим область нестабильности
        P, K = self._find_unstable_region(P, K, z, components)

        F = None
        trivial_counter = 0
        
        for j in range(self.max_iter):
            # 3. Состав паровой фазы
            y = z * K
            y_sum = np.sum(y)
            y_normalized = y / y_sum
            y_dict = {c: float(y_normalized[i]) for i, c in enumerate(components)}
            vapor_comp = self.composition.new_composition(y_dict)
            vapor_comp.T = self.T

            # 4. Расчет УРС
            eos_L = BrusilovskiyEOS(self.composition, P, self.T)
            eos_V = BrusilovskiyEOS(vapor_comp, P, self.T)

            eos_L.calc_eos()
            eos_V.calc_eos()

            Z_L = np.min(eos_L.real_roots_eos)
            Z_V = np.max(eos_V.real_roots_eos)

            # Проверка на тривиальное решение
            if np.isclose(Z_L, Z_V, rtol=1e-4):
                trivial_counter += 1
                logger.warning(f"Iter {j}: Trivial solution detected (Z_L ≈ Z_V = {Z_L:.6f})")
                
                if trivial_counter >= 3:
                    logger.info(f"Iter {j}: Escaping trivial solution via stability test")
                    P_new, K_new = self._escape_trivial_solution(P, K, z, components)
                    
                    if P_new < P * 0.95:  # Если давление значительно снизилось
                        P = P_new
                        if K_new is not None:
                            K = K_new
                        trivial_counter = 0
                        continue
                    else:
                        logger.warning("Cannot escape trivial solution. Stopping.")
                        break

            ln_phi_L = eos_L.get_fugacity_coef_vector_by_root(Z_L)
            ln_phi_V = eos_V.get_fugacity_coef_vector_by_root(Z_V)

            # 5. Новые K-факторы
            ln_K_new = ln_phi_L - ln_phi_V
            K_new = np.exp(ln_K_new)

            # 6. Функция F
            F = np.sum(z * K_new) - 1.0

            if abs(F) < self.tol:
                logger.info(f"Bubble point converged in {j+1} iterations. P = {P:.4f}, F = {F:.2e}")
                self.result_P = P
                self.result_F = F
                self.converged = True
                self.iterations_done = j + 1
                return P

            # 7. Производные
            eos_L._z_factor = Z_L
            eos_L._dz_dp = None
            eos_L._dlogphi_dp = None
            dlnphiL_dP = np.array([eos_L.calc_d_log_phi_i_dp(c) for c in components])

            eos_V._z_factor = Z_V
            eos_V._dz_dp = None
            eos_V._dlogphi_dp = None
            dlnphiV_dP = np.array([eos_V.calc_d_log_phi_i_dp(c) for c in components])

            dF_dP = np.sum(z * K_new * (dlnphiL_dP - dlnphiV_dP))

            if dF_dP > 0:
                logger.debug(f"Iter {j}: Sign correction for dF/dP. Original: {dF_dP:.2e} -> {-abs(dF_dP):.2e}")
                dF_dP = -abs(dF_dP)

            if abs(dF_dP) < 1e-15:
                logger.warning(f"dF/dP too small at iteration {j}. Stopping.")
                break

            # 8. Обновление давления
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

        logger.warning(f"Max iterations ({self.max_iter}) reached. P = {P:.4f}, F = {F:.2e}")
        self.result_P = P
        self.result_F = F
        self.converged = False
        self.iterations_done = self.max_iter
        return P



    def _find_unstable_region(self, P: float, K: np.ndarray, z: np.ndarray, components: list):
        """
        Проверяет стабильность при начальном P. Если система стабильна (P выше bubble point),
        использует бисекцию для нахождения давления, где система становится нестабильной.
        """
        stability_test = TwoPhaseStabilityTest(self.composition, P, self.T)
        stability_test.calculate_phase_stability()
        
        if not stability_test.stable:
            logger.info(f"Initial P={P:.4f} is unstable. Using stability test K-values.")
            if stability_test.k_values_for_flash is not None:
                K_new = np.array([stability_test.k_values_for_flash[c] for c in components])
                return P, K_new
            return P, K
        
        # Система стабильна → P слишком высокое. Ищем P, где система нестабильна.
        logger.info(f"Initial P={P:.4f} is stable (above bubble point). Searching for unstable region...")
        
        P_low = P * 0.1  # Нижняя граница (заведомо нестабильная)
        P_high = P       # Верхняя граница (стабильная)
        
        K_unstable = None
        
        # Бисекция для нахождения границы стабильности
        for i in range(30):
            P_mid = (P_low + P_high) / 2.0
            
            st = TwoPhaseStabilityTest(self.composition, P_mid, self.T)
            st.calculate_phase_stability()
            
            if st.stable:
                P_high = P_mid
            else:
                P_low = P_mid
                # Сохраняем K-факторы из нестабильной области
                if st.k_values_for_flash is not None:
                    K_unstable = np.array([st.k_values_for_flash[c] for c in components])
            
            if (P_high - P_low) / P_high < 0.05:  # 5% точность
                break
        
        P_new = P_low  # Берём давление из нестабильной области
        
        if K_unstable is not None:
            logger.info(f"Found unstable region at P={P_new:.4f}. Using stability test K-values.")
            return P_new, K_unstable
        else:
            logger.info(f"Found unstable region at P={P_new:.4f}. Using Wilson K-values.")
            # Пересчитываем K Вильсона при новом P
            comp_data = self.composition.composition_data
            Pc = np.array([comp_data['critical_pressure'][c] for c in components])
            Tc = np.array([comp_data['critical_temperature'][c] for c in components])
            omega = np.array([comp_data['acentric_factor'][c] for c in components])
            exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
            K_new = (Pc / P_new) * exp_term
            return P_new, K_new

    def _escape_trivial_solution(self, P: float, K: np.ndarray, z: np.ndarray, components: list):
        """
        Пытается выйти из тривиального решения, снижая давление и используя тест стабильности.
        """
        stability_test = TwoPhaseStabilityTest(self.composition, P, self.T)
        stability_test.calculate_phase_stability()
        
        if stability_test.stable:
            # Система стабильна → снижаем P
            P_new = P * 0.85
            logger.info(f"Stable at P={P:.4f}, reducing to P={P_new:.4f}")
            
            # Пересчитываем K Вильсона при новом P
            comp_data = self.composition.composition_data
            Pc = np.array([comp_data['critical_pressure'][c] for c in components])
            Tc = np.array([comp_data['critical_temperature'][c] for c in components])
            omega = np.array([comp_data['acentric_factor'][c] for c in components])
            exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
            K_new = (Pc / P_new) * exp_term
            
            return P_new, K_new
        else:
            # Система нестабильна → используем K из теста стабильности
            if stability_test.k_values_for_flash is not None:
                K_new = np.array([stability_test.k_values_for_flash[c] for c in components])
                logger.info(f"Unstable at P={P:.4f}. Using stability test K-values.")
                return P, K_new
            else:
                # Fallback: снижаем P
                P_new = P * 0.85
                comp_data = self.composition.composition_data
                Pc = np.array([comp_data['critical_pressure'][c] for c in components])
                Tc = np.array([comp_data['critical_temperature'][c] for c in components])
                omega = np.array([comp_data['acentric_factor'][c] for c in components])
                exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
                K_new = (Pc / P_new) * exp_term
                return P_new, K_new


# Для обратной совместимости
def calculate_bubble_point_pressure(
    composition: Composition,
    T: float,
    P_guess: float = None,
    max_iter: int = 200,
    tol: float = 1e-13,
) -> float:
    calculator = BubblePointCalculator(
        composition=composition,
        T=T,
        P_guess=P_guess,
        max_iter=max_iter,
        tol=tol,
    )
    return calculator.calculate()