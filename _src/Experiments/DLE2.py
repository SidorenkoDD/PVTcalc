from _src.Composition.CompositionV2 import Composition
from _src.Utils.Conditions import Conditions
from _src.VLE.Flash import Flash
from _src.PhaseDiagram.new_methodv2 import SaturationPressure
from _src.Utils.Conditions import Conditions, StandardConditions
from _src.Utils.Errors import LenthMissMatchError
import numpy as np


class DLE:
    def __init__(self, composition: Composition, pressure_arr_bar: list,
                 reservoir_pressure_bar: float, reservoir_temperature_c: float):
        self.composition = composition
        self.pressure_arr = pressure_arr_bar
        self.reservoir_pressure = reservoir_pressure_bar
        self.reservoir_temperature = reservoir_temperature_c

    @staticmethod
    def _is_descending(arr):
        if all(arr[i] < arr[i-1] for i in range(1, len(arr))) is False:
            arr.sort(reverse=True)
        else:
            pass


    @staticmethod
    def _check_length_arr(arr1, arr2):
        if len(arr1) == len(arr2):
            return True
        else:
            return False


    def _is_psat_in_pressure_arr(self):
        ...


    def _is_pres_in_pressure_arr(self):
        ...


    def _calculate_stage(self, p_bar, t_c):
        current_conditions = Conditions(p_bar, t_c)
        flash_object = Flash(self.composition, current_conditions)
        return flash_object.calculate()

    def _calculate_bo(self, liq_vol: np.ndarray, fl_arr: np.ndarray) -> np.ndarray:
        # Векторизованная замена None/NaN на 1.0
        fl_arr = np.where(np.isnan(fl_arr) | (fl_arr == None), 1.0, fl_arr).astype(float)
        liq_vol = np.array(liq_vol, dtype=float)
        
        # Накопительное произведение (cumulative product)
        # Сдвигаем на 1, чтобы первая итерация умножалась на 1, как в исходном коде
        cumulative_product = np.cumprod(np.concatenate(([1.0], fl_arr[:-1])))
        
        corrected_vol = liq_vol * cumulative_product
        self.oil_residual_volume = corrected_vol[-1]
        
        return corrected_vol / self.oil_residual_volume


    def _gas_vol_to_stc(self, p_stage: np.ndarray, t_stage: np.ndarray, 
                        z_stage: np.ndarray, v_stage: np.ndarray, z_stc: float) -> np.ndarray:
        return (p_stage * v_stage * z_stc * self.T_STC) / (self.P_STC * z_stage * t_stage)


    def _calculate_rs(self, p_arr: np.ndarray, z_arr: np.ndarray, t_arr: np.ndarray, 
                      gas_vol_arr: np.ndarray, fl_arr: np.ndarray, p_sat: float) -> np.ndarray:
        
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


    def calculate(self):
            result = []
            
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


            saturation_conditions = Conditions(saturation_pressure, self.reservoir_temperature)
            saturation_flash_object = Flash(self.composition, saturation_conditions)
            saturation_flash_result = saturation_flash_object.calculate()

            result.append(saturation_flash_result)
            self._liquid_molar_fractions = saturation_flash_result.liquid_composition

            ## ОСНОВНОЙ ЦИКЛ РАСЧЕТА ПО СТУПЕНЯМ
            for stage_pressure in self.pressure_arr:
                # Обновляем состав объекта композиции перед новым флешем
                self.composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)
                
                self._stage_result = self._calculate_stage(stage_pressure, self.reservoir_temperature)
                result.append(self._stage_result)
                
                # Снова безопасно берем состав жидкости для следующей итерации
                self._liquid_molar_fractions = self._stage_result.liquid_composition

            ## Расчет на стандартные условия
            stc_conditions = StandardConditions()
            # Передаем финальный состав жидкости на стандартные условия
            final_composition = self.composition.new_composition(self._liquid_molar_fractions, deep_copy=True)
            final_composition.T = stc_conditions.t
            
            stc_flash_object = Flash(final_composition, stc_conditions)
            stc_flash_result = stc_flash_object.calculate()
            result.append(stc_flash_result)

            return result
