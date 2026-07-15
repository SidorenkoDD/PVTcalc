"""
Дифференциация состава при постоянном составе (Constant Composition Expansion).

В отличие от `DLE` (где выделившийся газ удаляется), в CCE состав всей
системы остаётся постоянным на всех ступенях — считается общий (не
"остаточный жидкостный") объём. Используется для получения зависимости
объёма (`v/v_res`, `v/v_sat`) и сжимаемости от давления.
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.Utils.Conditions import Conditions
from calc_core.VLE.Flash import Flash
from calc_core.PhaseEnvelope.new_methodv2 import SaturationPressure
from calc_core.Utils.Conditions import Conditions, StandardConditions
from calc_core.Utils.Errors import LenthMissMatchError
import numpy as np
import pandas as pd
import numpy as np
from typing import List
from calc_core.VLE.Flash import FlashResult
from joblib import Parallel, delayed
from scipy.interpolate import UnivariateSpline

logger = logging.getLogger(__name__)


class CCE:
    """Один прогон CCE для заданного состава и сетки давлений при постоянной температуре."""

    def __init__(self, composition: Composition, pressure_arr_bar : list, reservoir_temperature:float):
        """
        Parameters
        ----------
        composition : Composition
            Состав (с уже посчитанными `composition_data`). Побочный эффект:
            `composition.T` сразу выставляется в `reservoir_temperature`
            (в K — параметр называется "температура", но подставляется как есть,
            без конвертации из °C — вызывающий код должен сам передать K).
        pressure_arr_bar : list
            Ступени давления, бар. Мутируется в `_calc_saturation_pressure`
            (в конец дописывается найденное P_sat).
        reservoir_temperature : float
            Температура эксперимента, K (см. caveat выше про `composition.T`).
        """
        self.composition = composition
        self.pressure_arr = pressure_arr_bar
        self.reservoir_temperature = reservoir_temperature
        self.composition.T = self.reservoir_temperature

    def _calc_stage(self, stage_pressure:float):
        """
        Флэш на одной ступени давления, на **копии** исходного состава
        (`deep_copy=True`) — в отличие от `calculate()`, здесь `self.composition`
        не переиспользуется и не мутируется. Используется в `calculate_parallel`.

        Parameters
        ----------
        stage_pressure : float
            Давление ступени, бар.

        Returns
        -------
        FlashResult
        """
        comp = self.composition.new_composition(self.composition.composition, deep_copy=True)
        condition_object = Conditions(stage_pressure, self.reservoir_temperature)
        flash_object = Flash(comp, condition_object)
        return flash_object.calculate()

    def _calc_saturation_pressure(self):
        """
        Находит давление насыщения (`PhaseEnvelope.new_methodv2.SaturationPressure`)
        для текущего состава при `reservoir_temperature + 273.15` — **обратите
        внимание**: здесь добавляется +273.15 поверх `reservoir_temperature`,
        то есть предполагается, что `reservoir_temperature` уже в Кельвинах
        (иначе результат будет вдвойне сконвертирован). Результат сохраняется
        в `self.saturation_pressure` и дописывается в `self.pressure_arr`.
        """
        sat_pressure_obj = SaturationPressure(self.composition, self.reservoir_temperature + 273.15)
        self.saturation_pressure = sat_pressure_obj.sp_convergence_loop()
        self.pressure_arr.append(self.saturation_pressure)


    def calculate(self):
        """
        Главная точка входа: последовательный (не параллельный, см.
        `calculate_parallel`) прогон CCE.

        Порядок: находим P_sat (дописывается в `self.pressure_arr`) → флэш на
        каждой ступени `self.pressure_arr` (на **одном и том же**
        `self.composition`, без копирования — если это важно для вас, см.
        `_calc_stage`, который копирует) → векторизация в DataFrame
        (`_vectorize_dle_results`) → `v/v_res`, `v/v_sat`, `Compressibility`
        как дополнительные колонки.

        Returns
        -------
        pd.DataFrame
            По одной строке на ступень давления (включая P_sat), отсортировано
            по убыванию давления. Колонки — свойства фаз (см.
            `_vectorize_dle_results`) + `v/v_res`, `v/v_sat`, `Compressibility`.
        """
        logger.info("CCE: старт, %d ступеней, T=%s°C", len(self.pressure_arr), self.reservoir_temperature)
        self._calc_saturation_pressure()
        logger.info("CCE: P_sat=%.4f бар", self.saturation_pressure)
        result = []


        # sat_pressure_obj = SaturationPressure(self.composition, self.reservoir_temperature + 273.15)
        # self.saturation_pressure = sat_pressure_obj.sp_convergence_loop()

        # saturation_conditions = Conditions(self.saturation_pressure, self.reservoir_temperature)
        # saturation_flash_object = Flash(self.composition, saturation_conditions)
        # saturation_flash_result = saturation_flash_object.calculate()
        # result.append(saturation_flash_result)

        for stage_num, stage_pressure in enumerate(self.pressure_arr, start=1):
            logger.debug("CCE: ступень %d/%d, P=%s бар", stage_num, len(self.pressure_arr), stage_pressure)
            stage_conditions = Conditions(stage_pressure, self.reservoir_temperature)
            stage_pressure_flash_object = Flash(self.composition, stage_conditions)
            result.append(stage_pressure_flash_object.calculate())

        df_res = self._vectorize_dle_results(result)
        self._calculate_v_d_vpres(df_res)
        self._calculate_v_d_vsat(df_res)
        self._calculate_compressibility(df_res)
        logger.info("CCE: завершено, %d точек в результате", len(result))
        return df_res
    

    def calculate_parallel(self):
        """
        То же, что `calculate()`, но ступени считаются параллельно через
        `joblib.Parallel` (`_calc_stage`, на независимых копиях состава).
        В отличие от `calculate()` — **не** делает векторизацию/постобработку
        (`v/v_res` и т.п.), возвращает сырой список `FlashResult`. По заметке
        в ноутбуке ("CCE в параллель — неэффективно") этот путь на практике
        оказался медленнее последовательного (накладные расходы на
        сериализацию `Composition` в отдельные процессы).

        Returns
        -------
        list[FlashResult]
        """
        self._calc_saturation_pressure()
        grid_res = Parallel(n_jobs=-1, backend="loky")(delayed(self._calc_stage)(P) for P in self.pressure_arr)
        return grid_res


    def _calculate_v_d_vpres(self, df):
        """
        Добавляет колонку `v/v_res` — отношение объёма к объёму на
        максимальном давлении (пластовые условия).

        **Известный баг**: условие `type(df['liquid_molar_volume']) != None`
        сравнивает тип Series с `None` — всегда `True`, поэтому ветка `else`
        (расчёт по `vapor_molar_volume`) фактически недостижима, даже если
        `liquid_molar_volume` весь `NaN` (однофазный газовый случай).

        Parameters
        ----------
        df : pd.DataFrame
            Модифицируется на месте.
        """
        if type(df['liquid_molar_volume']) != None:
            df['v/v_res'] = df['liquid_molar_volume'] / df[df['pressure'] == max(self.pressure_arr)]['liquid_molar_volume'].iloc[0]
        else:
            df['v/v_res'] = df['vapor_molar_volume'] / df[df['pressure'] == max(self.pressure_arr)]['vapor_molar_volume'].iloc[0]

    def _calculate_v_d_vsat(self, df):
        """
        Добавляет колонку `v/v_sat` — отношение объёма к объёму на давлении
        насыщения. Тот же баг с недостижимой веткой `else`, что и в
        `_calculate_v_d_vpres`.

        Parameters
        ----------
        df : pd.DataFrame
            Модифицируется на месте.
        """
        if type(df['liquid_molar_volume']) != None:
            df['v/v_sat'] = df['liquid_molar_volume'] / df[df['pressure'] == self.saturation_pressure]['liquid_molar_volume'].iloc[0]
        else:
            df['v/v_sat'] = df['vapor_molar_volume'] / df[df['pressure'] == self.saturation_pressure]['vapor_molar_volume'].iloc[0]



    def _calculate_compressibility(self, df, use_poly=False):
        """
        Добавляет колонку `Compressibility` — коэффициент сжимаемости
        `-(1/V)(dV/dP)`, посчитанный конечными разностями (`np.gradient` +
        минимальный из соседних объёмов в знаменателе) либо аналитически
        через полином 3-й степени, аппроксимирующий V(P) (`use_poly=True`).
        Для этого же есть отдельный публичный метод `calculate_compressibility`
        с другим API (принимает массивы напрямую, а не DataFrame состояния) —
        не связаны друг с другом, дублирующая функциональность.

        Parameters
        ----------
        df : pd.DataFrame
            Должен содержать `liquid_molar_volume`/`pressure`. Модифицируется на месте.
        use_poly : bool, optional
            `False` (по умолчанию) — конечно-разностный метод; `True` — через полином.

        Raises
        ------
        ValueError
            Если в `df` меньше 2 строк.
        """
        # 1. Извлекаем данные как NumPy массивы
        vol = df['liquid_molar_volume'].to_numpy(dtype=float)
        pres = df['pressure'].to_numpy(dtype=float)
        n = len(vol)
        
        # Защита от слишком коротких массивов (нужно минимум 2 точки для разности)
        if n < 2:
            raise ValueError("Для расчета сжимаемости конечно-разностным методом нужно минимум 2 точки.")
        
        # 2. Инициализируем массив
        compressibility = np.full(n, np.nan, dtype=float)
        
        # 3. Основной расчет
        if use_poly:
            # --- ВАРИАНТ С ПОЛИНОМОМ ---
            degree = 3
            coeffs = np.polyfit(pres, vol, degree)
            poly = np.poly1d(coeffs)
            poly_deriv = np.polyder(poly)
            
            # Для полинома понятие "минимум из двух точек" неприменимо, 
            # так как производная берется аналитически строго в точке pres[i].
            compressibility = np.abs((1.0 / vol) * poly_deriv(pres))
            
        else:
            # --- КЛАССИЧЕСКИЙ ВАРИАНТ С np.gradient и MIN(V1, V2) ---
            
            # Шаг 3.1: Считаем производную dV/dP во всех точках
            dV_dP = np.gradient(vol, pres)

            # Шаг 3.2: Формируем массив минимальных объемов для каждой схемы разности
            V_min = np.empty(n)
            
            # Для первой точки (прямая разность: используются индексы 0 и 1)
            V_min[0] = np.minimum(vol[0], vol[1])
            
            # Для внутренних точек (центральная разность: используются индексы i-1 и i+1)
            # vol[:-2] - это массив от 0 до n-3
            # vol[2:]  - это массив от 2 до n-1
            V_min[1:-1] = np.minimum(vol[:-2], vol[2:])
            
            # Для последней точки (обратная разность: используются индексы n-2 и n-1)
            V_min[-1] = np.minimum(vol[-2], vol[-1])
            
            # Шаг 3.3: Расчет сжимаемости
            # np.abs гарантирует положительный результат независимо от сортировки давления
            with np.errstate(divide='ignore', invalid='ignore'):
                compressibility = np.abs((1.0 / V_min) * dV_dP)

        df['Compressibility'] = compressibility


    def calculate_compressibility(self, df, pressures, volumes, method='average'):
        """
        Рассчитывает коэффициент сжимаемости между соседними точками.
        
        Parameters:
        pressures (list or np.array): Массив давлений (по убыванию).
        volumes (list or np.array): Согласованный массив объемов.
        method (str): Метод выбора объема в знаменателе.
                    'average' - среднее (V1 + V2) / 2 (Рекомендуется)
                    'initial' - начальное V1
                    'final'   - конечное V2
                    
        Returns:
        np.array: Массив значений сжимаемости (длина N-1).
        np.array: Массив средних давлений для каждой точки (для построения графиков).
        """
        P = np.asarray(pressures, dtype=float)
        V = np.asarray(volumes, dtype=float)
        
        if len(P) != len(V):
            raise ValueError("Массивы давлений и объемов должны иметь одинаковую длину.")
        if len(P) < 2:
            raise ValueError("Для расчета нужно минимум 2 точки.")
            
        # Рассчитываем разности (ΔP и ΔV)
        # np.diff вычисляет X[i+1] - X[i]
        dP = np.diff(P) 
        dV = np.diff(V)
        
        # Выбираем опорный объем (V_ref) в зависимости от метода
        if method == 'initial':
            V_ref = V[:-1]      # V1 (объемы от 0 до N-2)
        elif method == 'final':
            V_ref = V[1:]       # V2 (объемы от 1 до N-1)
        elif method == 'average':
            V_ref = (V[:-1] + V[1:]) / 2.0  # Среднее арифметическое
        else:
            raise ValueError("Неверный метод. Выберите 'average', 'initial' или 'final'.")
        
        # Защита от деления на ноль, если где-то давление не изменилось
        with np.errstate(divide='ignore', invalid='ignore'):
            beta = - (1.0 / V_ref) * (dV / dP)
            
        # Давления, соответствующие интервалам (середины интервалов для корректного графика)
        P_mid = (P[:-1] + P[1:]) / 2.0
        
        df['Compressibility'] = beta
        return beta, P_mid
        
    @staticmethod
    def _vectorize_dle_results(flash_results: List['FlashResult']) -> pd.DataFrame:
        """
        Преобразует список объектов FlashResult в векторизованный Pandas DataFrame.
        Скалярные свойства разворачиваются в колонки, составы остаются как dict.
        """
        records = []
        
        # 1. Динамически собираем все возможные ключи свойств из всех результатов
        # Это делает функцию устойчивой к изменениям в FluidPropertiesCalculator
        all_prop_keys = set()
        for res in flash_results:
            if res.vapor.properties:
                all_prop_keys.update(res.vapor.properties.keys())
            if res.liquid.properties:
                all_prop_keys.update(res.liquid.properties.keys())
                
        prop_keys = sorted(list(all_prop_keys)) 
        # Например: ['density', 'molar_density', 'molecular_weight', 'molar_volume', 'viscosity']

        # 2. Проходим по каждому шагу и формируем строку данных
        for res in flash_results:
            record = {
                # Базовые условия
                'pressure': res.pressure,
                'temperature': res.temperature,
                'is_two_phase': res.is_two_phase,
                
                # Мольные доли фаз (всегда есть, даже если 0.0)
                'vapor_mole_frac': res.vapor.mole_fraction,
                'liquid_mole_frac': res.liquid.mole_fraction,
                
                # Составы оставляем как есть (объекты dict), как вы и просили
                'vapor_composition': res.vapor.composition,
                'liquid_composition': res.liquid.composition,
            }
            
            # 3. Векторизуем свойства динамически
            for key in prop_keys:
                # Для газа: берем значение или np.nan, если словарь пуст или ключа нет
                v_prop = res.vapor.properties.get(key, np.nan) if res.vapor.properties else np.nan
                record[f'vapor_{key}'] = v_prop
                
                # Для жидкости: аналогично
                l_prop = res.liquid.properties.get(key, np.nan) if res.liquid.properties else np.nan
                record[f'liquid_{key}'] = l_prop
                
            records.append(record)
            
        # 4. Создаем DataFrame. Pandas автоматически определит типы:
        # float64 для чисел, bool для is_two_phase, object для словарей составов
        df = pd.DataFrame(records)
        
        # =====================================================================
        # ГАРАНТИЯ ПОРЯДКА: 
        # 1. Насильно присваиваем индекс строго по порядку следования в исходном массиве
        # df.index = range(len(df))
        # 2. Сортируем по этому индексу (для абсолютной гарантии) и сбрасываем 
        #    его в чистый числовой ряд (0, 1, 2...), убирая его из данных
        df = df.sort_values(by = 'pressure', ascending= False)
        # =====================================================================
        
        # Опционально: если в будущем потребуется сортировка по давлению, 
        # это нужно будет делать ПОСЛЕ гарантии исходного порядка, если это имеет смысл для вашей логики.
        # df = df.sort_values('pressure', ascending=False).reset_index(drop=True)
        
        return df