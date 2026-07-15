import numpy as np

from calc_core.Experiments.BaseExperiment import PVTExperiment
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.PhaseEnvelope.SaturationPressure_den_v import SaturationPointCalculator
from calc_core.PhaseEnvelope.CriticalProperties_den_v import CriticalPointCalculator
from calc_core.Utils.FluidPropertiesCalculator import FluidPropertiesCalculator
from calc_core.VLE.PhaseEquilibriumNewton import PhaseEquilibriumNewton
from calc_core.Utils.Constants import CONSTANT_R
from calc_core.Composition.Composition import Composition


class CVD(PVTExperiment):
    def __init__(self, composition: Composition):
        self._composition = composition
        self.result = {}

        self._components = tuple(composition.composition.keys())
        self._component_index = {c: i for i, c in enumerate(self._components)}
        self._nc = len(self._components)

        self._z0 = np.fromiter(
            (composition.composition[c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )

        self._molar_mass = np.fromiter(
            (composition.composition_data['molar_mass'][c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )

        self._c5_plus_flag = np.fromiter(
            (composition.composition_data['c5_plus_flag'][c] for c in self._components),
            dtype=bool,
            count=self._nc,
        )

    # =====================================================================================
    # ВСПОМОГАТЕЛЬНЫЕ МЕТОДЫ
    # =====================================================================================

    def _array_to_dict(self, arr: np.ndarray):
        return {c: float(arr[i]) for i, c in enumerate(self._components)}

    def _dict_to_array(self, dct: dict):
        return np.fromiter(
            (dct[c] for c in self._components),
            dtype=np.float64,
            count=self._nc,
        )

    def _eval_pc5_plus_array(self, composition_arr: np.ndarray) -> float:
        zi_Mi = np.sum(composition_arr[self._c5_plus_flag] * self._molar_mass[self._c5_plus_flag])
        return float(zi_Mi / 0.02404)

    def _find_critical_point(self):
        solver = CriticalPointCalculator(composition=self._composition)

        try:
            Tc, Pc = solver.calculate()
        except Exception as e:
            raise RuntimeError(f'Критическая точка для смеси не найдена: {e}')

        return Tc, Pc

    # =====================================================================================
    # ПОИСК Pdew
    # =====================================================================================

    def _find_dew_point(self, T: float, p_arr: np.ndarray):
        p_arr = np.sort(np.asarray(p_arr, dtype=np.float64))[::-1]

        p_upper = float(p_arr[0])
        p_lower = float(p_upper)

        stab_obj = TwoPhaseStabilityTest(self._composition, p_upper, T)
        stab_obj.calculate_phase_stability()

        if not stab_obj.stable:
            while not stab_obj.stable:
                p_lower = p_upper
                p_upper *= 1.1
                stab_obj = TwoPhaseStabilityTest(self._composition, p_upper, T)
                stab_obj.calculate_phase_stability()

            stab_obj = TwoPhaseStabilityTest(self._composition, p_lower, T)
            stab_obj.calculate_phase_stability()

        else:
            for p_candidate in p_arr[1:]:
                p_lower = float(p_candidate)

                stab_obj = TwoPhaseStabilityTest(self._composition, p_lower, T)
                stab_obj.calculate_phase_stability()

                if not stab_obj.stable:
                    break

                p_upper = p_lower

            if np.isclose(p_upper, p_lower):
                return None, None

        k_vals = stab_obj.k_vals_for_sat_point_calculation

        try:
            sat_point_obj = SaturationPointCalculator(self._composition)
            p_dew, iphase_comp = sat_point_obj.calculate(
                T, p_upper, p_lower, k_vals, sattype='dew')
            return float(p_dew), iphase_comp

        except Exception as e:
            raise RuntimeError(f'Давление начала ретроградной конденсации не найдено: {e}')

    # =====================================================================================
    # ОСНОВНОЙ РАСЧЕТ CVD
    # =====================================================================================

    def calculate(self, T: float, p_arr: list):
        self.result = {}

        p_arr = np.sort(np.asarray(p_arr, dtype=np.float64))
        p_dew, iphase_comp = self._find_dew_point(T, p_arr)

        if p_dew is not None:
            p_arr_full = np.unique(np.append(p_arr, p_dew))
            p_arr_full = np.sort(p_arr_full)[::-1]

            dew_idx = np.where(np.isclose(p_arr_full, p_dew, rtol=1e-12, atol=1e-12))[0]
            if dew_idx.size == 0:
                raise RuntimeError('Не удалось определить индекс давления насыщения в массиве давлений.')
            dew_idx = int(dew_idx[0])

            p_above_p_dew = p_arr_full[:dew_idx]
            p_below_p_dew = p_arr_full[dew_idx + 1:]
        else:
            p_above_p_dew = np.sort(p_arr)[::-1]
            p_below_p_dew = np.array([], dtype=np.float64)

        N_mole_mixture = 1.0
        mixture_composition_arr = self._z0.copy()
        mixture_composition_dict = self._array_to_dict(mixture_composition_arr)
        composition_properties = self._composition.composition_data

        # -----------------------------------------------------------------
        # Выше Pdew: однофазная область
        # -----------------------------------------------------------------
        for p in p_above_p_dew:
            res_p = {}
            mix_comp = self._composition.new_composition(mixture_composition_dict)

            gas_eos = BrusilovskiyEOS(composition=mix_comp, p=float(p), t=T)
            gas_eos.calc_eos()

            gas_properties = FluidPropertiesCalculator(
                mixture_composition_dict,
                composition_properties,
                gas_eos,
                float(p),
                T
            )

            res_p['gas_composition'] = mixture_composition_dict
            res_p['liquid_composition'] = None
            res_p['mixture_composition'] = mixture_composition_dict
            res_p['gas_viscosity'] = gas_properties.viscosity
            res_p['gas_z'] = gas_properties.z_shift
            res_p['gas_density'] = gas_properties.density
            res_p['gas_molar_volume'] = gas_properties.molar_volume
            res_p['gas_volume'] = N_mole_mixture * gas_properties.molar_volume
            res_p['liquid_viscosity'] = 0.0
            res_p['liquid_z'] = 0.0
            res_p['liquid_density'] = 0.0
            res_p['liquid_molar_volume'] = 0.0
            res_p['liquid_volume'] = 0.0
            res_p['Fl'] = 0.0
            res_p['Fv'] = 1.0
            res_p['pc5_plus'] = self._eval_pc5_plus_array(mixture_composition_arr)

            self.result[float(p)] = res_p

        # -----------------------------------------------------------------
        # Непосредственно в Pdew
        # -----------------------------------------------------------------
        if p_dew is not None:
            mix_comp = self._composition.new_composition(mixture_composition_dict)

            gas_eos = BrusilovskiyEOS(composition=mix_comp, p=float(p_dew), t=T)
            gas_eos.calc_eos()

            # iphase_comp_arr = self._dict_to_array(iphase_comp)
            liquid_composition = self._composition.new_composition(iphase_comp)
            liquid_eos = BrusilovskiyEOS(composition=liquid_composition, p=float(p_dew), t=T)
            liquid_eos.calc_eos()

            gas_properties = FluidPropertiesCalculator(
                mixture_composition_dict,
                composition_properties,
                gas_eos,
                float(p_dew),
                T
            )
            liquid_properties = FluidPropertiesCalculator(
                iphase_comp,
                composition_properties,
                liquid_eos,
                float(p_dew),
                T
            )

            V0 = N_mole_mixture * gas_properties.molar_volume

            self.result[float(p_dew)] = {}
            self.result[float(p_dew)]['gas_composition'] = mixture_composition_dict
            self.result[float(p_dew)]['liquid_composition'] = iphase_comp
            self.result[float(p_dew)]['mixture_composition'] = mixture_composition_dict
            self.result[float(p_dew)]['gas_viscosity'] = gas_properties.viscosity
            self.result[float(p_dew)]['gas_z'] = gas_properties.z_shift
            self.result[float(p_dew)]['gas_density'] = gas_properties.density
            self.result[float(p_dew)]['gas_molar_volume'] = gas_properties.molar_volume
            self.result[float(p_dew)]['gas_volume'] = V0
            self.result[float(p_dew)]['liquid_viscosity'] = liquid_properties.viscosity
            self.result[float(p_dew)]['liquid_z'] = liquid_properties.z_shift
            self.result[float(p_dew)]['liquid_density'] = liquid_properties.density
            self.result[float(p_dew)]['liquid_molar_volume'] = liquid_properties.molar_volume
            self.result[float(p_dew)]['liquid_volume'] = 0.0
            self.result[float(p_dew)]['Fl'] = 0.0
            self.result[float(p_dew)]['Fv'] = 1.0
            self.result[float(p_dew)]['pc5_plus'] = self._eval_pc5_plus_array(mixture_composition_arr)

            # -------------------------------------------------------------
            # Ниже Pdew: двухфазная область
            # -------------------------------------------------------------
            for p in p_below_p_dew:
                res_p = {}

                mix = self._composition.new_composition(mixture_composition_dict)

                stab_obj = TwoPhaseStabilityTest(mix, float(p), T)
                stab_obj.calculate_phase_stability()

                if stab_obj.stable:
                    continue

                phase_equilibrium = PhaseEquilibriumNewton(
                    mix,
                    float(p),
                    T,
                    stab_obj.k_values_for_flash
                )
                phase_equilibrium.find_solve_loop()

                liquid_composition_dict = phase_equilibrium.xi_l
                gas_composition_dict = phase_equilibrium.yi_v

                liquid_composition_arr = self._dict_to_array(liquid_composition_dict)
                gas_composition_arr = self._dict_to_array(gas_composition_dict)

                liquid_properties = FluidPropertiesCalculator(
                    liquid_composition_dict,
                    composition_properties,
                    phase_equilibrium.eos_liquid,
                    float(p),
                    T
                )
                gas_properties = FluidPropertiesCalculator(
                    gas_composition_dict,
                    composition_properties,
                    phase_equilibrium.eos_vapour,
                    float(p),
                    T
                )

                Fv = float(phase_equilibrium.fv)
                Fl = 1.0 - Fv
                Zv = gas_properties.z_shift

                Nl = N_mole_mixture * Fl
                Vl = Nl * liquid_properties.molar_volume
                Vv = V0 - Vl
                Nv = Vv * float(p) / (Zv * CONSTANT_R * T)

                N_mole_mixture = Nv + Nl

                res_p['gas_composition'] = gas_composition_dict
                res_p['liquid_composition'] = liquid_composition_dict
                res_p['mixture_composition'] = mixture_composition_dict
                res_p['gas_viscosity'] = gas_properties.viscosity
                res_p['gas_z'] = gas_properties.z_shift
                res_p['gas_density'] = gas_properties.density
                res_p['gas_molar_volume'] = gas_properties.molar_volume
                res_p['gas_volume'] = Vv
                res_p['liquid_viscosity'] = liquid_properties.viscosity
                res_p['liquid_z'] = liquid_properties.z_shift
                res_p['liquid_density'] = liquid_properties.density
                res_p['liquid_molar_volume'] = liquid_properties.molar_volume
                res_p['liquid_volume'] = Vl / V0 * 100.0
                res_p['Fl'] = Fl
                res_p['Fv'] = Fv
                res_p['pc5_plus'] = self._eval_pc5_plus_array(gas_composition_arr)

                self.result[float(p)] = res_p

                mixture_composition_arr = (
                    Nl * liquid_composition_arr + Nv * gas_composition_arr
                ) / N_mole_mixture
                mixture_composition_dict = self._array_to_dict(mixture_composition_arr)