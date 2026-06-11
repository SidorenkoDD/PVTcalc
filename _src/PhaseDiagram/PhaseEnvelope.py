import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Dict, Optional
import logging
from _src.Composition.CompositionV2 import Composition
from _src.PhaseDiagram.BubblePointPressure import BubblePointCalculator
from _src.PhaseDiagram.DewPressure import DewPointCalculator

logger = logging.getLogger(__name__)


class PhaseDiagramCalculator:
    """
    Класс для расчета фазовой диаграммы (P-T диаграммы) углеводородной смеси.
    """

    def __init__(
        self,
        composition: Composition,
        T_min: float = None,
        T_max: float = None,
        n_points: int = 20,
    ):
        self.composition = composition
        self.n_points = n_points
        
        comp_data = composition.composition_data
        z = np.array([composition.composition[c] for c in composition.composition.keys()])
        Tc_components = np.array([comp_data['critical_temperature'][c] for c in composition.composition.keys()])
        
        Tc_pseudo = np.sum(z * Tc_components)
        
        self.T_min = T_min if T_min is not None else 250.0
        self.T_max = T_max if T_max is not None else Tc_pseudo * 1.2
        self.Tc_pseudo = Tc_pseudo
        
        self.bubble_points = None
        self.dew_points_lower = None
        self.dew_points_upper = None
        self.critical_point = None
        
    def calculate_phase_envelope(
        self,
        T_range: np.ndarray = None,
        calculate_critical: bool = True,
        verbose: bool = True,
    ) -> Dict:
        """
        Расчет фазовой огибающей.
        """
        if T_range is None:
            T_range = np.linspace(self.T_min, self.T_max, self.n_points)
        
        bubble_T = []
        bubble_P = []
        dew_T_lower = []
        dew_P_lower = []
        dew_T_upper = []
        dew_P_upper = []
        
        logger.info(f"Starting phase envelope calculation for {len(T_range)} temperature points...")
        
        # Сначала находим критическую точку (приближенно)
        T_crit_approx = self.Tc_pseudo
        
        # Для передачи P_bubble как P_guess для upper dew point
        last_P_bubble = None
        
        for i, T in enumerate(T_range):
            if verbose:
                logger.info(f"Calculating point {i+1}/{len(T_range)} at T = {T:.2f} K ({T-273.15:.2f} °C)")
            
            self.composition.T = T
            
            # Bubble point (существует при T < T_crit)
            P_bubble = None
            try:
                bubble_calc = BubblePointCalculator(
                    composition=self.composition,
                    T=T,
                    max_iter=200,
                    tol=1e-12,
                )
                P_bubble = bubble_calc.calculate()
                
                if bubble_calc.converged and P_bubble is not None and P_bubble > 0:
                    bubble_T.append(T)
                    bubble_P.append(P_bubble)
                    last_P_bubble = P_bubble
                    if verbose:
                        logger.info(f"  Bubble point: P = {P_bubble:.4f} bar")
            except Exception as e:
                logger.warning(f"  Bubble point calculation failed at T={T:.2f} K: {e}")
            
            # Dew point: выбираем режим в зависимости от температуры
            if T < T_crit_approx:
                # При T < T_crit: считаем normal dew point (как upper)
                try:
                    dew_calc = DewPointCalculator(
                        composition=self.composition,
                        T=T,
                        dew_point_type='upper',
                        P_guess=last_P_bubble,  # Используем P_bubble как начальное приближение
                        max_iter=200,
                        tol=1e-12,
                    )
                    P_dew = dew_calc.calculate()
                    
                    if dew_calc.converged and P_dew is not None and P_dew > 0:
                        dew_T_upper.append(T)
                        dew_P_upper.append(P_dew)
                        if verbose:
                            logger.info(f"  Upper dew point (T < T_crit): P = {P_dew:.4f} bar")
                except Exception as e:
                    logger.warning(f"  Upper dew point calculation failed at T={T:.2f} K: {e}")
            else:
                # При T > T_crit: считаем upper и lower dew point
                # Upper dew point: используем P_bubble как P_guess
                try:
                    dew_calc_upper = DewPointCalculator(
                        composition=self.composition,
                        T=T,
                        dew_point_type='upper',
                        P_guess=last_P_bubble,  # Ключевое: передаем P_bubble как начальное приближение
                        max_iter=200,
                        tol=1e-12,
                    )
                    P_dew_upper = dew_calc_upper.calculate()
                    
                    if dew_calc_upper.converged and P_dew_upper is not None and P_dew_upper > 0:
                        dew_T_upper.append(T)
                        dew_P_upper.append(P_dew_upper)
                        if verbose:
                            logger.info(f"  Upper dew point (T > T_crit): P = {P_dew_upper:.4f} bar")
                except Exception as e:
                    logger.warning(f"  Upper dew point calculation failed at T={T:.2f} K: {e}")
                
                # Lower dew point: низкое начальное P
                try:
                    dew_calc_lower = DewPointCalculator(
                        composition=self.composition,
                        T=T,
                        dew_point_type='lower',
                        max_iter=200,
                        tol=1e-12,
                    )
                    P_dew_lower = dew_calc_lower.calculate()
                    
                    if dew_calc_lower.converged and P_dew_lower is not None and P_dew_lower > 0:
                        dew_T_lower.append(T)
                        dew_P_lower.append(P_dew_lower)
                        if verbose:
                            logger.info(f"  Lower dew point (T > T_crit): P = {P_dew_lower:.4f} bar")
                except Exception as e:
                    logger.warning(f"  Lower dew point calculation failed at T={T:.2f} K: {e}")
        
        self.bubble_points = {
            'T': np.array(bubble_T),
            'P': np.array(bubble_P),
            'T_C': np.array(bubble_T) - 273.15,
        }
        
        self.dew_points_lower = {
            'T': np.array(dew_T_lower),
            'P': np.array(dew_P_lower),
            'T_C': np.array(dew_T_lower) - 273.15,
        }
        
        self.dew_points_upper = {
            'T': np.array(dew_T_upper),
            'P': np.array(dew_P_upper),
            'T_C': np.array(dew_T_upper) - 273.15,
        }
        
        if calculate_critical:
            self.critical_point = self._find_critical_point()
        
        return {
            'bubble_point': self.bubble_points,
            'dew_point_lower': self.dew_points_lower,
            'dew_point_upper': self.dew_points_upper,
            'critical_point': self.critical_point,
        }
    
    def _find_critical_point(self) -> Dict:
        """
        Поиск критической точки.
        """
        if self.bubble_points is None or self.dew_points_upper is None:
            logger.warning("Phase envelope not calculated yet. Cannot find critical point.")
            return None
        
        P_bubble = self.bubble_points['P']
        P_dew_upper = self.dew_points_upper['P']
        T_bubble = self.bubble_points['T']
        T_dew_upper = self.dew_points_upper['T']
        
        if len(P_bubble) == 0 or len(P_dew_upper) == 0:
            logger.warning("Insufficient data to find critical point.")
            return None
        
        T_common = np.linspace(max(T_bubble.min(), T_dew_upper.min()), 
                               min(T_bubble.max(), T_dew_upper.max()), 100)
        
        P_bubble_interp = np.interp(T_common, T_bubble, P_bubble)
        P_dew_interp = np.interp(T_common, T_dew_upper, P_dew_upper)
        
        delta_P = np.abs(P_bubble_interp - P_dew_interp)
        
        min_idx = np.argmin(delta_P)
        
        T_crit = T_common[min_idx]
        P_crit = (P_bubble_interp[min_idx] + P_dew_interp[min_idx]) / 2.0
        
        logger.info(f"Critical point found: T = {T_crit:.2f} K ({T_crit-273.15:.2f} °C), P = {P_crit:.2f} bar")
        
        self.critical_point = {
            'T': T_crit,
            'P': P_crit,
            'T_C': T_crit - 273.15,
        }
        
        return self.critical_point
    
    def plot_phase_diagram(
        self,
        show: bool = True,
        save_path: str = None,
        figsize: Tuple[int, int] = (10, 8),
        title: str = None,
    ) -> plt.Figure:
        """
        Отрисовка фазовой диаграммы P-T.
        """
        if self.bubble_points is None:
            raise RuntimeError("Phase envelope not calculated. Call calculate_phase_envelope() first.")
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Bubble point curve
        if len(self.bubble_points['T']) > 0:
            ax.plot(
                self.bubble_points['T_C'], 
                self.bubble_points['P'], 
                'b-', 
                linewidth=2, 
                label='Bubble Point',
                marker='o',
                markersize=4,
            )
        
        # Upper dew point curve
        if len(self.dew_points_upper['T']) > 0:
            ax.plot(
                self.dew_points_upper['T_C'], 
                self.dew_points_upper['P'], 
                'r-', 
                linewidth=2, 
                label='Upper Dew Point',
                marker='^',
                markersize=4,
            )
        
        # Lower dew point curve
        if len(self.dew_points_lower['T']) > 0:
            ax.plot(
                self.dew_points_lower['T_C'], 
                self.dew_points_lower['P'], 
                'r--', 
                linewidth=2, 
                label='Lower Dew Point',
                marker='v',
                markersize=4,
            )
        
        # Critical point
        if self.critical_point is not None:
            ax.plot(
                self.critical_point['T_C'], 
                self.critical_point['P'], 
                'g*', 
                markersize=15, 
                label=f'Critical Point\n({self.critical_point["T_C"]:.1f}°C, {self.critical_point["P"]:.1f} bar)',
                zorder=5,
            )
        
        ax.set_xlabel('Temperature, °C', fontsize=12)
        ax.set_ylabel('Pressure, bar', fontsize=12)
        
        if title is None:
            title = 'Phase Diagram (P-T)'
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Заштриховываем двухфазную область
        self._shade_two_phase_region(ax)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Phase diagram saved to {save_path}")
        
        if show:
            plt.show()
        
        return fig
    
    def _shade_two_phase_region(self, ax):
        """
        Заштриховка двухфазной области.
        """
        # Область между bubble point и upper dew point (при T < T_crit)
        if len(self.bubble_points['T']) > 0 and len(self.dew_points_upper['T']) > 0:
            T_min = max(self.bubble_points['T_C'].min(), self.dew_points_upper['T_C'].min())
            T_max = min(self.bubble_points['T_C'].max(), self.dew_points_upper['T_C'].max())
            
            if T_max > T_min:
                T_fill = np.linspace(T_min, T_max, 100)
                P_bubble_fill = np.interp(T_fill, self.bubble_points['T_C'], self.bubble_points['P'])
                P_dew_fill = np.interp(T_fill, self.dew_points_upper['T_C'], self.dew_points_upper['P'])
                
                ax.fill_between(
                    T_fill, 
                    P_bubble_fill, 
                    P_dew_fill, 
                    alpha=0.2, 
                    color='gray',
                    label='Two-Phase Region (T < T_crit)',
                )
        
        # Область между lower и upper dew point (при T > T_crit, ретроградная конденсация)
        if len(self.dew_points_lower['T']) > 0 and len(self.dew_points_upper['T']) > 0:
            T_min = max(self.dew_points_lower['T_C'].min(), self.dew_points_upper['T_C'].min())
            T_max = min(self.dew_points_lower['T_C'].max(), self.dew_points_upper['T_C'].max())
            
            if T_max > T_min:
                T_fill = np.linspace(T_min, T_max, 100)
                P_dew_lower_fill = np.interp(T_fill, self.dew_points_lower['T_C'], self.dew_points_lower['P'])
                P_dew_upper_fill = np.interp(T_fill, self.dew_points_upper['T_C'], self.dew_points_upper['P'])
                
                ax.fill_between(
                    T_fill, 
                    P_dew_lower_fill, 
                    P_dew_upper_fill, 
                    alpha=0.2, 
                    color='orange',
                    label='Retrograde Condensation Region (T > T_crit)',
                )
    
    def get_phase_envelope_data(self) -> Dict:
        """
        Возвращает данные фазовой огибающей.
        """
        return {
            'bubble_point': self.bubble_points,
            'dew_point_lower': self.dew_points_lower,
            'dew_point_upper': self.dew_points_upper,
            'critical_point': self.critical_point,
            'composition': self.composition.composition,
            'T_range': {
                'min': self.T_min,
                'max': self.T_max,
                'min_C': self.T_min - 273.15,
                'max_C': self.T_max - 273.15,
            },
        }
