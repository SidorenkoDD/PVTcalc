"""
Многоступенчатая сепарация (промысловый сепараторный тест).

Структурно почти идентичен `DLE` (см. его module docstring) — то же
"пластовое условие → P_sat → цикл по ступеням → стандартные условия"
с переносом состава жидкости между ступенями, — но здесь у каждой ступени
своя **пара** (P, T), а не только давление при общей температуре.
Соответствующие методы (`_calculate_bo`, `_gas_vol_to_stc`, `_calculate_rs`,
`_vectorize_dle_results`) — почти дословные копии из `DLE.py`, см. пометку
в CLAUDE.md про дублирование кода между этими двумя файлами.
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.Utils.Conditions import Conditions
from calc_core.VLE.Flash import Flash
from calc_core.PhaseEnvelope.new_methodv2 import SaturationPressure
from calc_core.Utils.Conditions import Conditions, StandardConditions
from calc_core.Utils.Errors import LengthMismatchError
import numpy as np
import pandas as pd
import numpy as np
from typing import List
from calc_core.VLE.Flash import FlashResult

logger = logging.getLogger(__name__)

class SeparatorTest:
    """Один прогон многоступенчатой сепарации для заданного состава и сетки (P, T) по ступеням."""

    def __init__(self, composition : Composition,
                 pressure_arr_bar : list,
                 temperature_arr_c : list,
                 reservoir_pressure_bar : float,
                 reservoir_temperature_c : float):
        """
        Parameters
        ----------
        composition : Composition
            Пластовый состав.
        pressure_arr_bar : list
            Давления ступеней, бар.
        temperature_arr_c : list
            Температуры ступеней, °C — **должен быть той же длины**, что
            `pressure_arr_bar` (проверяется `_check_length_arr` в `calculate()`).
        reservoir_pressure_bar : float
            Пластовое давление, бар.
        reservoir_temperature_c : float
            Пластовая температура, °C.
        """


        self.composition = composition
        self.pressure_arr = pressure_arr_bar
        self.temperature_arr = temperature_arr_c
        self.reservoir_pressure = reservoir_pressure_bar
        self.reservoir_temperature = reservoir_temperature_c
        self.stc_conditions = StandardConditions()
    
    @staticmethod
    def _is_descending(arr):
        """Сортирует `arr` по убыванию на месте, если ещё не отсортирован. Нигде в классе не вызывается (мёртвый код), см. аналог в `DLE._is_descending`."""
        if all(arr[i] < arr[i-1] for i in range(1, len(arr))) is False:
            arr.sort(reverse=True)
        else:
            pass


    @staticmethod
    def _check_length_arr(arr1, arr2):
        """
        Проверяет, что два массива одинаковой длины (в `calculate()` —
        `pressure_arr` и `temperature_arr`). В отличие от `DLE._check_length_arr`
        (возвращает bool), здесь при несовпадении **бросает исключение**.

        Raises
        ------
        LenthMissMatchError
            Если длины не совпадают (опечатка в имени класса исключения,
            см. CLAUDE.md Known Issues).
        """
        if len(arr1) == len(arr2):
            return True
        else:
            raise LengthMismatchError


    def _is_psat_in_pressure_arr(self):
        """Не реализовано (заготовка, тело — `...`)."""
        ...


    def _is_pres_in_pressure_arr(self):
        """Не реализовано (заготовка, тело — `...`)."""
        ...


    def _calculate_stage(self, p_bar, t_c):
        """
        Флэш на текущем `self.composition` для одной ступени.

        Parameters
        ----------
        p_bar : float
        t_c : float
            Температура, °C.

        Returns
        -------
        FlashResult
        """
        current_conditions = Conditions(p_bar, t_c)
        flash_object = Flash(self.composition, current_conditions)
        return flash_object.calculate()


    def _calculate_bo(self, liq_vol: np.ndarray, fl_arr: np.ndarray) -> np.ndarray:
        """
        Объёмный коэффициент нефти Bo по ступеням. Логика идентична
        `DLE._calculate_bo` (см. её docstring за подробностями).

        Parameters
        ----------
        liq_vol : np.ndarray
        fl_arr : np.ndarray

        Returns
        -------
        np.ndarray
        """
        # Векторизованная замена None/NaN на 1.0
        fl_arr = np.where(np.isnan(fl_arr) | (fl_arr == None), 1.0, fl_arr).astype(float)
        liq_vol = np.array(liq_vol, dtype=float)
        liq_vol = liq_vol * fl_arr
        # Накопительное произведение (cumulative product)
        # Сдвигаем на 1, чтобы первая итерация умножалась на 1, как в исходном коде
        cumulative_product = np.cumprod(np.concatenate(([1.0], fl_arr[:-1])))
        
        corrected_vol = liq_vol * cumulative_product
        self.oil_residual_volume = corrected_vol[-1]
        
        return corrected_vol / self.oil_residual_volume


    def _gas_vol_to_stc(self, p_stage: np.ndarray, t_stage: np.ndarray,
                        z_stage: np.ndarray, v_stage: np.ndarray, z_stc: float) -> np.ndarray:
        """Пересчёт объёма газа со ступени на стандартные условия. Идентично `DLE._gas_vol_to_stc`."""
        return (p_stage * v_stage * z_stc * self.stc_conditions.t) / (self.stc_conditions.p * z_stage * t_stage)


    def _calculate_rs(self, p_arr: np.ndarray, z_arr: np.ndarray, t_arr: np.ndarray,
                      gas_vol_arr: np.ndarray, fl_arr: np.ndarray, p_sat: float) -> np.ndarray:
        """
        Газосодержание Rs по ступеням. Логика идентична `DLE._calculate_rs`
        (см. её docstring за подробностями); требует, чтобы `_calculate_bo`
        уже был вызван (использует `self.oil_residual_volume`).

        Parameters
        ----------
        p_arr, z_arr, t_arr, gas_vol_arr, fl_arr : np.ndarray
        p_sat : float

        Returns
        -------
        np.ndarray
        """

        fl_arr = np.where(np.isnan(fl_arr) | (fl_arr == None), 1.0, fl_arr).astype(float)
        gas_vol_arr = np.array(gas_vol_arr, dtype=float)
        p_arr = np.array(p_arr, dtype=float)
        
        # Шаг 1: Коррекция объема газа
        cumulative_product = np.cumprod(np.concatenate(([1.0], fl_arr[:-1])))
        corrected_vol_fl = gas_vol_arr * cumulative_product

        # Шаг 2: Конвертация в STC (векторизованно!)
        z_stc = float(z_arr[-1]) # Явное приведение к скаляру
        gas_vol_stc_arr = self._gas_vol_to_stc(p_arr, t_arr, z_arr, corrected_vol_fl, z_stc)

        # Шаг 3: Обнуление выше p_sat и на последней ступени
        gas_vol_stc_arr = np.where(p_arr >= p_sat, 0.0, gas_vol_stc_arr)
        gas_vol_stc_arr[-1] = 0.0

        # Шаг 4 и 5: Кумулятивная сумма и инверсия (векторизованно)
        cumulative_sum = np.cumsum(gas_vol_stc_arr)
        gas_stc_acc_reverted = cumulative_sum[-1] - cumulative_sum
        
        return gas_stc_acc_reverted / self.oil_residual_volume


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
        df.index = range(len(df))
        # 2. Сортируем по этому индексу (для абсолютной гарантии) и сбрасываем 
        #    его в чистый числовой ряд (0, 1, 2...), убирая его из данных
        df = df.sort_index().reset_index(drop=True)
        # =====================================================================
        
        # Опционально: если в будущем потребуется сортировка по давлению, 
        # это нужно будет делать ПОСЛЕ гарантии исходного порядка, если это имеет смысл для вашей логики.
        # df = df.sort_values('pressure', ascending=False).reset_index(drop=True)
        
        return df


    def calculate(self):
            """
            Главная точка входа: полный прогон сепараторного теста.

            Порядок идентичен `DLE.calculate()`, но цикл идёт по парам
            `zip(pressure_arr, temperature_arr)` вместо одного `pressure_arr`
            (и температура состава явно переставляется на каждой ступени:
            `self.composition.T = stage_temperature`).

            Returns
            -------
            list[FlashResult]
                Пластовая точка, P_sat, каждая ступень, стандартные условия —
                итого `len(pressure_arr) + 3`. Таблица `Bo`/`Rs` — в `self._dle_df`.

            Raises
            ------
            LenthMissMatchError
                Если `pressure_arr` и `temperature_arr` разной длины.
            """
            logger.info("SeparatorTest: старт, %d ступеней", len(self.pressure_arr))
            result = []
            self._check_length_arr(self.pressure_arr, self.temperature_arr)

            ## Расчет для пластовых условий
            reservoir_conditions = Conditions(self.reservoir_pressure, self.reservoir_temperature)
            reservoir_flash_object = Flash(self.composition, reservoir_conditions)
            reservoir_flash_result = reservoir_flash_object.calculate()
            result.append(reservoir_flash_result)
            self._liquid_molar_fractions = reservoir_flash_result.liquid_composition

            ## Расчет для Рнас

            self.composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)
            sat_pressure_obj = SaturationPressure(self.composition, reservoir_conditions.t)
            saturation_pressure = sat_pressure_obj.sp_convergence_loop()
            logger.info("SeparatorTest: P_sat=%.4f бар", saturation_pressure)

            saturation_conditions = Conditions(saturation_pressure, self.reservoir_temperature)
            saturation_flash_object = Flash(self.composition, saturation_conditions)
            saturation_flash_result = saturation_flash_object.calculate()

            result.append(saturation_flash_result)
            self._liquid_molar_fractions = saturation_flash_result.liquid_composition

            ## ОСНОВНОЙ ЦИКЛ РАСЧЕТА ПО СТУПЕНЯМ
            for stage_num, (stage_pressure, stage_temperature) in enumerate(zip(self.pressure_arr, self.temperature_arr), start=1):
                logger.debug(
                    "SeparatorTest: ступень %d/%d, P=%s бар, T=%s°C",
                    stage_num, len(self.pressure_arr), stage_pressure, stage_temperature,
                )
                # Обновляем состав объекта композиции перед новым флешем
                self.composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)
                self.composition.T = stage_temperature
                self._stage_result = self._calculate_stage(stage_pressure, stage_temperature)
                result.append(self._stage_result)

                # Снова безопасно берем состав жидкости для следующей итерации
                self._liquid_molar_fractions = self._stage_result.liquid_composition

            ## Расчет на стандартные условия
            # Передаем финальный состав жидкости на стандартные условия
            final_composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)
            final_composition.T = self.stc_conditions.t
            
            stc_flash_object = Flash(final_composition, self.stc_conditions)
            stc_flash_result = stc_flash_object.calculate()
            result.append(stc_flash_result)

            self._dle_df = self._vectorize_dle_results(result)

            self._bo = self._calculate_bo(self._dle_df['liquid_molar_volume'].to_numpy(), self._dle_df['liquid_mole_frac'].to_numpy())
            self._rs = self._calculate_rs(p_arr = self._dle_df['pressure'].to_numpy(),
                                          z_arr = self._dle_df['vapor_z'].to_numpy(), 
                                          t_arr = self._dle_df['temperature'].to_numpy(),
                                          gas_vol_arr = self._dle_df['vapor_molar_volume'].to_numpy() * self._dle_df['vapor_mole_frac'].to_numpy(),
                                          fl_arr = self._dle_df['liquid_mole_frac'].to_numpy(),
                                          p_sat = saturation_conditions.p)
            self._dle_df['Bo'] = self._bo
            self._dle_df['Rs'] = self._rs

            logger.info("SeparatorTest: завершено, %d точек в результате", len(result))
            return result
