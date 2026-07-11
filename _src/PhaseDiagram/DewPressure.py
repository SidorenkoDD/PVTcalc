import numpy as np
import logging
from _src.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from _src.Composition.Composition import Composition

logger = logging.getLogger(__name__)

class DewPointCalculator:
    """
    Класс для расчета давления начала конденсации (Dew Point Pressure)
    по алгоритму из Table 6.2 (Pedersen).
    
    dew_point_type:
      'upper' - верхнее давление начала конденсации (высокое P)
      'lower' - нижнее давление начала конденсации (низкое P)
    """

    def __init__(
        self,
        composition: Composition,
        T: float,
        dew_point_type: str = 'upper',
        P_guess: float = None,
        max_iter: int = 100,
        tol: float = 1e-10,
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

    def calculate(self) -> float:
        """
        Расчет давления начала конденсации по алгоритму Table 6.2.
        """
        logger.info(f"Starting dew point calculation at T={self.T:.2f} K ({self.T-273.15:.2f} °C), type={self.dew_point_type}")
        
        # Получаем данные компонент
        comp_data = self.composition.composition_data
        components = list(self.composition.composition.keys())
        z = np.array([self.composition.composition[c] for c in components])
        
        # Критические параметры для начального приближения
        Pc = np.array([comp_data['critical_pressure'][c] for c in components])
        Tc = np.array([comp_data['critical_temperature'][c] for c in components])
        omega = np.array([comp_data['acentric_factor'][c] for c in components])
        
        # Шаг 1: Начальное приближение давления
        if self.P_guess is not None:
            P = self.P_guess
            logger.info(f"Using provided initial pressure guess: P={P:.4f} bar")
        else:
            # Оцениваем начальное давление по Вильсону
            exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
            if self.dew_point_type == 'upper':
                P = np.sum(z * Pc * exp_term) * 1.5
            else:
                K_wilson = (Pc / 1.0) * exp_term
                P = 1.0 / np.sum(z / K_wilson)
            logger.info(f"Initial pressure estimate (Wilson): P={P:.4f} bar")
        
        # Шаг 2: Оцениваем K-факторы по уравнению Вильсона
        exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - Tc / self.T))
        K = (Pc / P) * exp_term
        logger.debug(f"Initial K-factors (Wilson): {dict(zip(components, K))}")
        
        # Основной цикл итераций (шаги 3-9)
        for j in range(self.max_iter):
            logger.debug(f"\n=== Iteration {j+1}/{self.max_iter} ===")
            logger.debug(f"Current pressure: P={P:.6f} bar")
            
            # Шаг 3: Оцениваем состав жидкой фазы: x_i = z_i / K_i
            x = z / K
            x_sum = np.sum(x)
            logger.debug(f"Liquid phase composition sum: sum(x)={x_sum:.6f}")
            
            if x_sum <= 0:
                logger.error(f"Non-positive liquid phase sum at iteration {j+1}. Terminating.")
                break
                
            # Нормализуем состав жидкой фазы
            x_normalized = x / x_sum
            x_dict = {c: float(x_normalized[i]) for i, c in enumerate(components)}
            
            # Создаем объект для жидкой фазы
            liquid_comp = self.composition.new_composition(x_dict, deep_copy=True)
            liquid_comp.T = self.T
            logger.debug(f"Liquid composition: {x_dict}")
            
            # Шаг 4: Рассчитываем коэффициенты летучести для паровой и жидкой фаз
            # Паровая фаза (исходный состав = feed composition)
            eos_V = BrusilovskiyEOS(self.composition, P, self.T)
            eos_V.calc_eos()
            
            # Жидкая фаза (рассчитанный состав)
            eos_L = BrusilovskiyEOS(liquid_comp, P, self.T)
            eos_L.calc_eos()
            
            # Выбираем корни: для пара - максимальный, для жидкости - минимальный
            Z_V = np.max(eos_V.real_roots_eos)
            Z_L = np.min(eos_L.real_roots_eos)
            
            logger.debug(f"Vapor Z-factor: {Z_V:.6f}")
            logger.debug(f"Liquid Z-factor: {Z_L:.6f}")
            logger.debug(f"Z-factor difference: {abs(Z_V - Z_L):.6e}")
            
            # Получаем логарифмы коэффициентов летучести
            ln_phi_V = eos_V.get_fugacity_coef_vector_by_root(Z_V)
            ln_phi_L = eos_L.get_fugacity_coef_vector_by_root(Z_L)
            
            logger.debug(f"ln(phi_V): {dict(zip(components, ln_phi_V))}")
            logger.debug(f"ln(phi_L): {dict(zip(components, ln_phi_L))}")
            
            # Шаг 5: Рассчитываем новые K-факторы: K_i = phi_i^L / phi_i^V
            ln_K_new = ln_phi_L - ln_phi_V
            K_new = np.exp(ln_K_new)
            
            logger.debug(f"New K-factors from fugacity: {dict(zip(components, K_new))}")
            
            # Шаг 6: Вычисляем целевую функцию F = sum(z_i / K_i) - 1
            F = np.sum(z / K_new) - 1.0
            
            logger.info(f"Iteration {j+1}: P={P:.6f} bar, F={F:.6e}")
            logger.debug(f"sum(z/K) = {np.sum(z / K_new):.10f}")
            
            # Проверка сходимости
            if abs(F) < self.tol:
                logger.info(f"Dew point converged in {j+1} iterations!")
                logger.info(f"Final pressure: P={P:.6f} bar")
                logger.info(f"Final F={F:.2e} (tol={self.tol:.2e})")
                
                self.result_P = P
                self.result_F = F
                self.converged = True
                self.iterations_done = j + 1
                return P
            
            # Шаг 7: Вычисляем производную dF/dP
            # dF/dP = sum( (z_i / K_i) * (d(ln_phi_i^V)/dP - d(ln_phi_i^L)/dP) )
            
            # Устанавливаем Z-факторы для расчета производных
            eos_V._z_factor = Z_V
            eos_L._z_factor = Z_L
            
            # Рассчитываем производные ln(phi_i) по давлению
            dlnphiV_dP = np.array([eos_V.calc_d_log_phi_i_dp(c) for c in components])
            dlnphiL_dP = np.array([eos_L.calc_d_log_phi_i_dp(c) for c in components])
            
            logger.debug(f"d(ln_phi_V)/dP: {dict(zip(components, dlnphiV_dP))}")
            logger.debug(f"d(ln_phi_L)/dP: {dict(zip(components, dlnphiL_dP))}")
            
            # Вычисляем dF/dP
            dF_dP = np.sum((z / K_new) * (dlnphiV_dP - dlnphiL_dP))
            
            logger.debug(f"dF/dP = {dF_dP:.6e}")
            
            # Проверка на малую производную
            if abs(dF_dP) < 1e-15:
                logger.warning(f"dF/dP too small ({dF_dP:.2e}) at iteration {j+1}. Terminating.")
                break
            
            # Шаг 8: Обновляем давление по методу Ньютона
            # P_new = P - F / (dF/dP)
            delta_P = -F / dF_dP
            
            logger.debug(f"Newton step: delta_P = -F/(dF/dP) = {-F:.6e} / {dF_dP:.6e} = {delta_P:.6e}")
            
            # Ограничиваем шаг (не более 50% от текущего давления)
            max_step = 0.5 * P
            if abs(delta_P) > max_step:
                delta_P = max_step * np.sign(delta_P)
                logger.debug(f"Step limited to {max_step:.6f} bar")
            
            P_new = P + delta_P
            
            # Проверка на отрицательное давление
            if P_new <= 0:
                logger.warning(f"New pressure would be negative ({P_new:.6f}). Reducing step.")
                P_new = P * 0.5
                delta_P = P_new - P
            
            logger.debug(f"Updated pressure: P_new = {P:.6f} + {delta_P:.6f} = {P_new:.6f} bar")
            
            # Обновляем давление и K-факторы
            P = P_new
            K = K_new
            
            # Шаг 9: Если не сошлось, возвращаемся к шагу 3
            if (j + 1) % 10 == 0:
                logger.info(f"Iteration {j+1}: P={P:.6f} bar, F={F:.6e}, dF/dP={dF_dP:.6e}")
        
        # Достигнуто максимальное число итераций
        logger.warning(f"Maximum iterations ({self.max_iter}) reached.")
        logger.warning(f"Final pressure: P={P:.6f} bar, F={F:.6e}")
        
        self.result_P = None
        self.result_F = F
        self.converged = False
        self.iterations_done = self.max_iter
        return None


# Функция для обратной совместимости
def calculate_dew_point_pressure(
    composition: Composition,
    T: float,
    dew_point_type: str = 'upper',
    P_guess: float = None,
    max_iter: int = 100,
    tol: float = 1e-10,
) -> float:
    """
    Calculate dew point pressure at given temperature.
    """
    calculator = DewPointCalculator(
        composition=composition,
        T=T,
        dew_point_type=dew_point_type,
        P_guess=P_guess,
        max_iter=max_iter,
        tol=tol,
    )
    return calculator.calculate()