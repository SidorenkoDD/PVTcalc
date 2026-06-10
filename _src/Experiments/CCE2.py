from _src.Composition.CompositionV2 import Composition
from _src.Utils.Conditions import Conditions
from _src.VLE.Flash import Flash
from _src.PhaseDiagram.new_methodv2 import SaturationPressure
from _src.Utils.Conditions import Conditions, StandardConditions
from _src.Utils.Errors import LenthMissMatchError
import numpy as np
import pandas as pd
import numpy as np
from typing import List
from _src.VLE.Flash import FlashResult
from joblib import Parallel, delayed




class CCE:

    def __init__(self, composition: Composition, pressure_arr_bar : list, reservoir_temperature:float):
        self.composition = composition
        self.pressure_arr = pressure_arr_bar
        self.reservoir_temperature = reservoir_temperature
        self.composition.T = self.reservoir_temperature

    def _calc_stage(self, stage_pressure:float):
        comp = self.composition.new_composition(self.composition.composition, deep_copy=True)
        condition_object = Conditions(stage_pressure, self.reservoir_temperature)
        flash_object = Flash(comp, condition_object)
        return flash_object.calculate()
    


    def calculate(self):
        result = []


        sat_pressure_obj = SaturationPressure(self.composition, self.reservoir_temperature + 273.15)
        self.saturation_pressure = sat_pressure_obj.sp_convergence_loop()

        saturation_conditions = Conditions(self.saturation_pressure, self.reservoir_temperature)
        saturation_flash_object = Flash(self.composition, saturation_conditions)
        saturation_flash_result = saturation_flash_object.calculate()

        result.append(saturation_flash_result)

        for stage_pressure in self.pressure_arr:
            stage_conditions = Conditions(stage_pressure, self.reservoir_temperature)
            stage_pressure_flash_object = Flash(self.composition, stage_conditions)
            result.append(stage_pressure_flash_object.calculate())
        
        df_res = self._vectorize_dle_results(result)
        self._calculate_v_d_vpres(df_res)
        self._calculate_v_d_vsat(df_res)
        return df_res
    

    def calculate_parallel(self):
        grid_res = Parallel(n_jobs=-1, backend="loky")(delayed(self._calc_stage(P))(P) for P in self.pressure_arr)
        return grid_res
    

    def _calculate_v_d_vpres(self, df):
        if type(df['liquid_molar_volume']) != None:
            df['v/v_res'] = df['liquid_molar_volume'] / df[df['pressure'] == max(self.pressure_arr)]['liquid_molar_volume'].iloc[0]
        else:
            df['v/v_res'] = df['vapor_molar_volume'] / df[df['pressure'] == max(self.pressure_arr)]['vapor_molar_volume'].iloc[0]

    def _calculate_v_d_vsat(self, df):
        if type(df['liquid_molar_volume']) != None:
            df['v/v_sat'] = df['liquid_molar_volume'] / df[df['pressure'] == self.saturation_pressure]['liquid_molar_volume'].iloc[0]
        else:
            df['v/v_sat'] = df['vapor_molar_volume'] / df[df['pressure'] == self.saturation_pressure]['vapor_molar_volume'].iloc[0]


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