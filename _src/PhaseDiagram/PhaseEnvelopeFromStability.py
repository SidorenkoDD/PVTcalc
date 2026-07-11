import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from _src.Composition.Composition import Composition
from _src.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest

class PhaseEnvelopeFromStability:  # Исправлена опечатка в названии
    def __init__(self, composition: Composition,
                 max_pressure: float = 800,
                 max_temperature: float = 800,
                 pressure_points: int = 30,
                 temperature_points: int = 30):
        self.composition = composition
        self.max_pressure = max_pressure
        self.max_temperature = max_temperature
        self.pressure_points = pressure_points
        self.temperature_points = temperature_points

        # Хранилища результатов
        self.result_pressure_arr = []
        self.result_temperature_arr = []
        self.result_stability_flag_arr = []

        # Сетки параметров
        self.pressure_array = np.linspace(1, max_pressure, pressure_points)
        self.temperature_array = np.linspace(273.15, max_temperature + 273.15, temperature_points)

    def _calc_single_point(self, p: float, t: float):
        """
        Чистая функция расчёта для одной точки. 
        Не изменяет self, а возвращает кортеж результатов.
        """
        composition = self.composition.new_composition(self.composition.composition, deep_copy=True)
        composition.T = t
        # Передаём t напрямую, чтобы избежать гонки за self.composition.T
        stab_test_obj = TwoPhaseStabilityTest(composition=composition, p=p, t=t)
        stab_test_obj.calculate_phase_stability()

        flag = 0.0 if stab_test_obj.stable else 1.0
        return p, t - 273.15, flag  # P, T_C, Flag

    def run_parallel(self, n_jobs: int = -1, backend: str = 'loky'):
        """
        Запускает расчёт по всей сетке P-T в параллельном режиме.
        """
        # 1. Создаём все комбинации (декартово произведение)
        P_grid, T_grid = np.meshgrid(self.pressure_array, self.temperature_array)
        points = list(zip(P_grid.ravel(), T_grid.ravel()))

        # 2. Параллельный запуск
        results = Parallel(n_jobs=n_jobs, backend=backend)(
            delayed(self._calc_single_point)(p, t) for p, t in points
        )

        # 3. Безопасная сборка результатов в атрибуты класса
        if results:
            self.result_pressure_arr = [r[0] for r in results]
            self.result_temperature_arr = [r[1] for r in results]
            self.result_stability_flag_arr = [r[2] for r in results]
        else:
            self.result_pressure_arr = []
            self.result_temperature_arr = []
            self.result_stability_flag_arr = []

    def plot(self):
        plt.scatter(x=self.result_temperature_arr, y=self.result_pressure_arr, c=self.result_stability_flag_arr)
        plt.xlabel('T, C')
        plt.ylabel('P, bar')
        plt.title('Фазовая диаграмма только через тест стабильности')
        plt.show()

