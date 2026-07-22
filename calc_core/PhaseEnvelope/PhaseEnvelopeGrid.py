import threading

import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

from calc_core.Composition.Composition import Composition
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.Utils.Cancellation import CancellationToken, ProgressCallback, report_progress


class PhaseEnvelopeGrid:
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

    def _calc_single_point(self, p: float, t: float, cancellation_token=None):
        """
        Чистая функция расчёта для одной точки. 
        Не изменяет self, а возвращает кортеж результатов.
        """
        composition = self.composition.new_composition(self.composition.composition, deep_copy=True)
        composition.T = t
        # Передаём t напрямую, чтобы избежать гонки за self.composition.T
        stab_test_obj = TwoPhaseStabilityTest(
            composition=composition, p=p, t=t,
            cancellation_token=cancellation_token,
        )
        stab_test_obj.calculate_phase_stability()

        flag = 0.0 if stab_test_obj.stable else 1.0
        return p, t - 273.15, flag  # P, T_C, Flag

    def run_parallel(
        self,
        n_jobs: int = -1,
        backend: str = 'loky',
        cancellation_token: CancellationToken | None = None,
        progress_callback: ProgressCallback | None = None,
    ):
        """
        Запускает расчёт по всей сетке P-T в параллельном режиме.
        """
        # 1. Создаём все комбинации (декартово произведение)
        P_grid, T_grid = np.meshgrid(self.pressure_array, self.temperature_array)
        points = list(zip(P_grid.ravel(), T_grid.ravel()))

        # Потоковый backend нужен только для отменяемого GUI-пути: обычный
        # `loky` запускает отдельные процессы и не видит threading.Event из
        # CancellationToken. Численный алгоритм и порядок результатов те же.
        if cancellation_token is None:
            results = Parallel(n_jobs=n_jobs, backend=backend)(
                delayed(self._calc_single_point)(p, t) for p, t in points
            )
        else:
            cancellation_token.throw_if_cancelled()
            completed = 0
            counter_lock = threading.Lock()

            def calculate_point(p, t):
                nonlocal completed
                cancellation_token.throw_if_cancelled()
                result = self._calc_single_point(p, t, cancellation_token)
                with counter_lock:
                    completed += 1
                    report_progress(
                        progress_callback, completed / len(points),
                        f"Grid point {completed}/{len(points)}",
                    )
                return result

            results = Parallel(n_jobs=n_jobs, backend='threading')(
                delayed(calculate_point)(p, t) for p, t in points
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

    def run(
        self,
        cancellation_token: CancellationToken,
        progress_callback: ProgressCallback | None = None,
    ):
        """Последовательный отменяемый путь для GUI."""
        P_grid, T_grid = np.meshgrid(self.pressure_array, self.temperature_array)
        points = list(zip(P_grid.ravel(), T_grid.ravel()))
        results = []
        total = max(1, len(points))
        for index, (p, t) in enumerate(points, start=1):
            cancellation_token.throw_if_cancelled()
            results.append(self._calc_single_point(p, t, cancellation_token))
            report_progress(progress_callback, index / total, f"Grid point {index}/{total}")
        self.result_pressure_arr = [r[0] for r in results]
        self.result_temperature_arr = [r[1] for r in results]
        self.result_stability_flag_arr = [r[2] for r in results]

    def plot(self):
        plt.scatter(x=self.result_temperature_arr, y=self.result_pressure_arr, c=self.result_stability_flag_arr)
        plt.xlabel('T, C')
        plt.ylabel('P, bar')
        plt.title('Фазовая диаграмма только через тест стабильности')
        plt.show()

