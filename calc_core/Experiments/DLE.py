"""
Дифференциальная конденсация (Differential Liberation Experiment).

Стандартный PVT-эксперимент: пластовые условия → находим давление насыщения
(`PhaseEnvelope.new_methodv2.SaturationPressure`) → ступенчато снижаем
давление при постоянной температуре, на каждой ступени выделившийся газ
считается полностью удалённым из системы (в отличие от CCE), в расчёт на
следующей ступени идёт только жидкость (`composition.new_composition(...)`
с составом жидкости предыдущей ступени) → замер на стандартных условиях.
Итог — таблица `Bo`/`Rs` по ступеням давления (`self._dle_df`).
"""

import logging

from calc_core.Composition.Composition import Composition
from calc_core.Utils.Cancellation import CancellationToken, ProgressCallback, report_progress
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

class DLE:
    """Один прогон дифференциальной конденсации для заданного состава и сетки давлений."""

    def __init__(self, composition: Composition, pressure_arr_bar: list,
                 reservoir_pressure_bar: float, reservoir_temperature_c: float):
        """
        Parameters
        ----------
        composition : Composition
            Пластовый состав (с уже посчитанными `composition_data`).
        pressure_arr_bar : list
            Ступени давления, бар (порядок — как передан; на строго
            убывающий порядок код не проверяет, несмотря на наличие
            неиспользуемого хелпера `_is_descending`).
        reservoir_pressure_bar : float
            Пластовое давление, бар — первая точка расчёта.
        reservoir_temperature_c : float
            Пластовая (и всех ступеней — DLE изотермический) температура, °C.
        """
        self.composition = composition
        self.pressure_arr = pressure_arr_bar
        self.reservoir_pressure = reservoir_pressure_bar
        self.reservoir_temperature = reservoir_temperature_c
        self.stc_conditions = StandardConditions()


    
    def _is_descending(arr):
        """
        Сортирует `arr` по убыванию на месте, если он ещё не отсортирован.
        Не статический метод, но объявлен без `self` — вызов как `self._is_descending(...)`
        передал бы `self` в `arr`. Фактически нигде в классе не вызывается (мёртвый код).
        """
        if all(arr[i] < arr[i-1] for i in range(1, len(arr))) is False:
            arr.sort(reverse=True)
        else:
            pass


    @staticmethod
    def _check_length_arr(arr1, arr2):
        """Проверяет, что два массива одинаковой длины. Возвращает bool (в отличие от `SeparatorTest._check_length_arr`, который бросает исключение)."""
        if len(arr1) == len(arr2):
            return True
        else:
            return False


    def _is_psat_in_pressure_arr(self):
        """Не реализовано (заготовка, тело — `...`)."""
        ...


    def _is_pres_in_pressure_arr(self):
        """Не реализовано (заготовка, тело — `...`)."""
        ...


    def _calculate_stage(self, p_bar, t_c, cancellation_token=None, progress_callback=None):
        """
        Флэш на текущем `self.composition` при заданных P/T одной ступени.

        Parameters
        ----------
        p_bar : float
        t_c : float
            Температура, °C (конвертируется в K внутри `Conditions`).

        Returns
        -------
        FlashResult
        """
        current_conditions = Conditions(p_bar, t_c)
        flash_object = Flash(
            self.composition, current_conditions,
            cancellation_token=cancellation_token,
            progress_callback=progress_callback,
        )
        return flash_object.calculate()


    def _calculate_bo(self, liq_vol: np.ndarray, fl_arr: np.ndarray) -> np.ndarray:
        """
        Объёмный коэффициент нефти Bo по ступеням — отношение объёма жидкости
        на ступени (с накопленной поправкой на долю жидкости, оставшейся
        после каждой предыдущей ступени) к объёму на последней (стандартной)
        ступени (`self.oil_residual_volume`).

        Parameters
        ----------
        liq_vol : np.ndarray
            Молярный объём жидкости по ступеням (`liquid_molar_volume` из `self._dle_df`).
        fl_arr : np.ndarray
            Мольная доля жидкости по ступеням (`liquid_mole_frac`); `NaN`
            трактуется как 1.0 (вся ступень — жидкость).

        Returns
        -------
        np.ndarray
            Bo по ступеням, тот же порядок, что вход. Заодно сохраняет
            `self.oil_residual_volume` (объём на последней ступени) —
            используется потом и в `_calculate_rs`.
        """
        # 1. Безопасная замена None/NaN: сначала приводим к float (None станет np.nan), затем заменяем nan на 1.0
        fl_arr = np.array(fl_arr, dtype=float)
        fl_arr = np.where(np.isnan(fl_arr), 1.0, fl_arr)
        
        liq_vol = np.array(liq_vol, dtype=float)
        liq_vol = liq_vol * fl_arr
        # 2. Накопительное произведение со сдвигом (полный аналог цикла)
        cumulative_product = np.cumprod(np.concatenate(([1.0], fl_arr[:-1])))
        
        corrected_vol = liq_vol * cumulative_product
        self.oil_residual_volume = corrected_vol[-1]
        
        return corrected_vol / self.oil_residual_volume


    def _gas_vol_to_stc(self, p_stage: np.ndarray, t_stage: np.ndarray,
                        z_stage: np.ndarray, v_stage: np.ndarray, z_stc: float) -> np.ndarray:
        """
        Пересчёт объёма газа со ступени на стандартные условия (STC) через
        закон реального газа: `V_stc = P·V·Z_stc·T_stc / (P_stc·Z·T)`.

        Parameters
        ----------
        p_stage, t_stage, z_stage, v_stage : np.ndarray
            Давление, температура, Z-фактор и объём газа на ступени.
        z_stc : float
            Z-фактор газа на стандартных условиях (обычно с последней ступени).

        Returns
        -------
        np.ndarray
        """
        # Убедитесь, что self.stc_conditions.t == 293.14 (или 293.15) и p == 0.101325
        return (p_stage * v_stage * z_stc * self.stc_conditions.t) / \
               (self.stc_conditions.p * z_stage * t_stage)


    def _calculate_rs(self, p_arr: np.ndarray, z_arr: np.ndarray, t_arr: np.ndarray,
                      gas_vol_arr: np.ndarray, fl_arr: np.ndarray, p_sat: float) -> np.ndarray:
        """
        Газосодержание Rs по ступеням — накопленный (от текущей ступени до
        конца) объём выделившегося газа в STC, делённый на
        `self.oil_residual_volume`. Ступени выше `p_sat` (однофазные, газ ещё
        не выделился) обнуляются перед накоплением.

        Requires
        --------
        Вызывается после `_calculate_bo` (использует `self.oil_residual_volume`).

        Parameters
        ----------
        p_arr, z_arr, t_arr : np.ndarray
            Давление, Z-фактор пара, температура по ступеням.
        gas_vol_arr : np.ndarray
            Молярный объём выделившегося газа по ступеням (до пересчёта в STC).
        fl_arr : np.ndarray
            Мольная доля жидкости по ступеням (для накопленной поправки, как в `_calculate_bo`).
        p_sat : float
            Давление насыщения — ступени с `P >= p_sat` считаются безгазовыми.

        Returns
        -------
        np.ndarray
            Rs по ступеням, тот же порядок, что вход.
        """

        # 1. Безопасная обработка fl_arr (как в _calculate_bo)
        fl_arr = np.array(fl_arr, dtype=float)
        fl_arr = np.where(np.isnan(fl_arr), 1.0, fl_arr)
        
        # gas_vol_arr приводим к float, но НЕ заменяем NaN на 1.0, чтобы сохранить логику оригинала
        gas_vol_arr = np.array(gas_vol_arr, dtype=float)
        p_arr = np.array(p_arr, dtype=float)
        
        # 2. Коррекция объема газа (векторизованный cumprod)
        cumulative_product = np.cumprod(np.concatenate(([1.0], fl_arr[:-1])))
        corrected_vol_fl = gas_vol_arr * cumulative_product

        # 3. Конвертация в STC (векторизованно)
        z_stc = float(z_arr[-1]) 
        gas_vol_stc_arr = self._gas_vol_to_stc(p_arr, t_arr, z_arr, corrected_vol_fl, z_stc)

        # 4. Обнуление выше p_sat и на последней ступени
        gas_vol_stc_arr = np.where(p_arr >= p_sat, 0.0, gas_vol_stc_arr)
        # gas_vol_stc_arr[-1] = 0.0

        # 5. Кумулятивная сумма и инверсия (векторизованно)
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


    def calculate(self, cancellation_token: CancellationToken | None = None,
                  progress_callback: ProgressCallback | None = None):
            """
            Главная точка входа: полный прогон дифференциальной конденсации.

            Порядок: флэш на пластовых условиях → находим P_sat
            (`PhaseEnvelope.new_methodv2.SaturationPressure`) → флэш на P_sat →
            цикл по `self.pressure_arr` (флэш на каждой ступени, состав
            жидкости переходит на следующую через `new_composition(...,
            deep_copy=True)`) → флэш на стандартных условиях → векторизация
            в `self._dle_df` (`_vectorize_dle_results`) → расчёт `Bo`/`Rs`
            (`_calculate_bo`/`_calculate_rs`), добавленных туда же колонками.

            Returns
            -------
            list[FlashResult]
                По одному на каждую точку: пластовая, P_sat, каждая ступень
                `pressure_arr`, стандартные условия — итого `len(pressure_arr) + 3`.
                Таблица `Bo`/`Rs` доступна отдельно в `self._dle_df`.
            """
            logger.info("DLE: старт, %d ступеней, T=%s°C", len(self.pressure_arr), self.reservoir_temperature)
            if cancellation_token is not None:
                cancellation_token.throw_if_cancelled()
            result = []
            self.composition.T = self.reservoir_temperature
            ## Расчет для пластовых условий
            reservoir_conditions = Conditions(self.reservoir_pressure, self.reservoir_temperature)
            reservoir_flash_object = Flash(
                self.composition, reservoir_conditions,
                cancellation_token=cancellation_token,
                progress_callback=progress_callback,
            )
            reservoir_flash_result = reservoir_flash_object.calculate()
            result.append(reservoir_flash_result)
            self._liquid_molar_fractions = reservoir_flash_result.liquid_composition

            ## Расчет для Рнас

            self.composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)
            sat_pressure_obj = SaturationPressure(self.composition, reservoir_conditions.t)
            saturation_pressure = sat_pressure_obj.sp_convergence_loop()
            logger.info("DLE: P_sat=%.4f бар", saturation_pressure)

            saturation_conditions = Conditions(saturation_pressure, self.reservoir_temperature)
            saturation_flash_object = Flash(
                self.composition, saturation_conditions,
                cancellation_token=cancellation_token,
                progress_callback=progress_callback,
            )
            saturation_flash_result = saturation_flash_object.calculate()

            result.append(saturation_flash_result)
            self._liquid_molar_fractions = saturation_flash_result.liquid_composition

            ## ОСНОВНОЙ ЦИКЛ РАСЧЕТА ПО СТУПЕНЯМ
            for stage_num, stage_pressure in enumerate(self.pressure_arr, start=1):
                if cancellation_token is not None:
                    cancellation_token.throw_if_cancelled()
                report_progress(
                    progress_callback,
                    0.2 + 0.65 * (stage_num - 1) / max(1, len(self.pressure_arr)),
                    f"DLE stage {stage_num}/{len(self.pressure_arr)}",
                )
                logger.debug("DLE: ступень %d/%d, P=%s бар", stage_num, len(self.pressure_arr), stage_pressure)
                # Обновляем состав объекта композиции перед новым флешем
                self.composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)

                self._stage_result = self._calculate_stage(
                    stage_pressure, self.reservoir_temperature,
                    cancellation_token, progress_callback,
                )
                result.append(self._stage_result)

                # Снова безопасно берем состав жидкости для следующей итерации
                self._liquid_molar_fractions = self._stage_result.liquid_composition

            ## Расчет на стандартные условия
            # Передаем финальный состав жидкости на стандартные условия
            final_composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)
            final_composition.T = self.stc_conditions.t
            
            stc_flash_object = Flash(
                final_composition, self.stc_conditions,
                cancellation_token=cancellation_token,
                progress_callback=progress_callback,
            )
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

            logger.info("DLE: завершено, %d точек в результате", len(result))
            return result
