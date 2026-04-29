import math

# from .BasePhaseStability import PhaseStabilityTest
from _src.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from _src.Composition.CompositionV2 import Composition
# from ...calculations.EOS.EOSFactory import EOSFactory
from _src.Utils.Constants import TOL_TWO_PHASE_STABILITY_CONVERGENCE, TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION
from _src.Utils.Errors import StopIterationError


# class TwoPhaseStabilityTest(PhaseStabilityTest):
class TwoPhaseStabilityTest:

    def __init__(self, composition: Composition, p: float, t: float):

        self.p = p
        self.t = t

        self.stable = None
        self.k_values_liquid = None
        self.k_values_vapour = None
        self.k_values_for_flash = None
        self.S_l = None
        self.S_v = None
        self.convergence_l = False
        self.convergence_v = False
        self.convergence_trivial_solution_l = False
        self.convergence_trivial_solution_v = False

        # self.eos = EOSFactory.create_eos(eos_name)
        self._composition = composition
        self.zi = composition.composition
        self.composition_data = composition.composition_data

    @staticmethod
    def _calc_k_values_wilson(p, T, composition: dict, composition_data: dict):
        k_values = {}

        for component in list(composition.keys()):
            acf = composition_data['acentric_factor'][component]
            pc = composition_data['critical_pressure'][component]
            Tc = composition_data['critical_temperature'][component]

            k_values[component] = math.exp(5.37 * (1 + acf) * (1 - (Tc / T))) / (p / pc)

        return k_values

    # Расчет мольных долей в газовой фазе
    @staticmethod
    def _calc_Yi_v(zi: dict, k_values_vapour: dict):
        Yi_v = {}
        for component in list(k_values_vapour.keys()):
            Yi_v[component] = zi[component] * k_values_vapour[component]
        return Yi_v

    # Расчет мольных долей в жидкой фазе
    @staticmethod
    def _calc_Xi_l(zi: dict, k_values_liquid: dict):
        Xi_l = {}
        for component in list(k_values_liquid.keys()):
            Xi_l[component] = zi[component] / k_values_liquid[component]

        return Xi_l

    # Расчет суммы мольных долей в газовой фазе
    @staticmethod
    def _calc_S_v(Yi_v: dict):
        return sum(list(Yi_v.values()))

    # Расчет суммы мольных долей в жидкой фазе
    @staticmethod
    def _calc_S_l(Xi_l: dict):
        return sum(list(Xi_l.values()))

    # Нормируем мольные доли газовой фазы
    @staticmethod
    def _normalize_mole_fractions_vapour(Yi_v: dict, S_v: float):
        normalized_mole_fractions_vapour = {}
        for component in list(Yi_v.keys()):
            normalized_mole_fractions_vapour[component] = Yi_v[component] / S_v
        return normalized_mole_fractions_vapour

    # Нормируем мольные доли для жидкой фазы
    @staticmethod
    def _normalize_mole_fractions_liquid(Xi_l: dict, S_l: float):
        normalized_mole_fractions_liquid = {}
        for component in list(Xi_l.keys()):
            normalized_mole_fractions_liquid[component] = Xi_l[component] / S_l
        return normalized_mole_fractions_liquid

    # Рассчитываем Ri для газовой фазы
    @staticmethod
    def _calc_ri_vapour(vapour_fugacities: dict, mixture_fugacities: dict, S_v: float):
        ri_vapour = {}
        for component in mixture_fugacities.keys():
            ri = math.exp(mixture_fugacities[component]) / (math.exp(vapour_fugacities[component]) * S_v)
            ri_vapour[component] = ri
        return ri_vapour

    # Рассчитываем Ri для жидкой фазы
    @staticmethod
    def _calc_ri_liquid(liquid_figacities: dict, mixture_fugacities: dict, S_l: float):
        ri_liquid = {}
        for component in mixture_fugacities.keys():
            ri = math.exp(liquid_figacities[component]) / math.exp(mixture_fugacities[component]) * S_l
            ri_liquid[component] = ri
        return ri_liquid

    @staticmethod
    def _update_k_values(k_values_old: dict, ri: dict):
        new_k_i = {}
        for component in k_values_old.keys():
            new_k_i[component] = k_values_old[component] * ri[component]

        return new_k_i

    @staticmethod
    def _check_convergence(ri: dict) -> bool:
        e = TOL_TWO_PHASE_STABILITY_CONVERGENCE
        sq_ri_to_sum = []

        for val in list(ri.values()):
            sq_ri_to_sum.append(math.pow(val - 1, 2))

        sum_sq_ri = sum(sq_ri_to_sum)

        if sum_sq_ri < e:
            return True

        else:
            return False

    @staticmethod
    def _check_trivial_solution(k_values: dict):
        e = TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION

        sq_log_ki_to_sum = []
        for ki_v in list(k_values.values()):
            sq_log_ki_to_sum.append(math.pow(math.log(ki_v), 2))

        if sum(sq_log_ki_to_sum) < e:
            return True

        return False

    def _interpetate_stability_analysis(self):

        if ((self.convergence_trivial_solution_v and self.convergence_trivial_solution_l) or
                ((self.S_v <= 1) and self.convergence_trivial_solution_l) or
                (self.convergence_trivial_solution_v and (self.S_l <= 1)) or
                ((self.S_v <= 1) and (self.S_l <= 1))):
            self.stable = True

        elif (self.S_v > 1) and self.convergence_trivial_solution_l:
            self.stable = False
            self.k_values_for_flash = self.k_values_vapour

        elif self.convergence_trivial_solution_v and (self.S_l > 1):
            self.stable = False
            self.k_values_for_flash = self.k_values_liquid

        elif (self.S_v > 1) and (self.S_l > 1):
            self.stable = False
            self.k_values_for_flash = {}

            for component in self.k_values_liquid.keys():
                k_v = self.k_values_vapour[component]
                k_l = self.k_values_liquid[component]
                self.k_values_for_flash[component] = k_v * k_l

        elif (self.S_v > 1) and (self.S_l <= 1):
            self.stable = False
            self.k_values_for_flash = self.k_values_vapour

        else:
            self.stable = False
            self.k_values_for_flash = self.k_values_liquid

    def calculate_phase_stability(self):
        p = self.p
        T = self.t

        # initial_eos = self.eos(composition=self._composition, p=p, t=T)
        initial_eos = BrusilovskiyEOS(composition=self._composition, p=p, t=T)
        initial_eos.calc_eos()
        mixture_fugacities = initial_eos.fugacities

        k_values_vapour_init = self._calc_k_values_wilson(p, T, self.zi, self.composition_data)
        k_values_liquid_init = self._calc_k_values_wilson(p, T, self.zi, self.composition_data)

        self._loop_vapour(p, T, k_values_vapour_init, mixture_fugacities)
        self._loop_liquid(p, T, k_values_liquid_init, mixture_fugacities)

        self._interpetate_stability_analysis()

    def _loop_vapour(self, p: float, T: float, k_values: dict, mixture_fugacities: dict):
        i = 0
        self.k_values_vapour = k_values

        while True:
            Yi_v = self._calc_Yi_v(self.zi, self.k_values_vapour)
            self.S_v = self._calc_S_v(Yi_v)
            self.yi_v = self._normalize_mole_fractions_vapour(Yi_v, self.S_v)

            new_composition = self._composition.new_composition(self.yi_v)
            # new_composition.T = T
            self.vapour_eos = BrusilovskiyEOS(composition=new_composition, p=p, t=T)
            # self.vapour_eos = self.eos(composition=self._composition.new_composition(self.yi_v), p=p, t=T)
            self.vapour_eos.calc_eos()
            vapour_fugacities = self.vapour_eos.fugacities

            ri_v = self._calc_ri_vapour(vapour_fugacities, mixture_fugacities, self.S_v)

            if self._check_convergence(ri_v):
                self.convergence_v = True
                break

            self.k_values_vapour = self._update_k_values(self.k_values_vapour, ri_v)

            if self._check_trivial_solution(self.k_values_vapour):
                self.convergence_trivial_solution_v = True
                break

            i += 1

            if i > 100000:
                raise StopIterationError('Число итераций теста стабильности превысило 100000')

    def _loop_liquid(self, p: float, T: float, k_values: dict, mixture_fugacities: dict):
        i = 0
        self.k_values_liquid = k_values

        while True:
            Xi_l = self._calc_Xi_l(self.zi, self.k_values_liquid)
            self.S_l = self._calc_S_l(Xi_l)
            self.xi_l = self._normalize_mole_fractions_vapour(Xi_l, self.S_l)

            new_composition = self._composition.new_composition(self.xi_l)
            # new_composition.T = T
            self.liquid_eos = BrusilovskiyEOS(composition=new_composition, p=p, t=T)
            # self.liquid_eos = self.eos(composition=self._composition.new_composition(self.xi_l), p=p, t=T)
            self.liquid_eos.calc_eos()
            liquid_fugacities = self.liquid_eos.fugacities

            ri_l = self._calc_ri_liquid(liquid_fugacities, mixture_fugacities, self.S_l)

            if self._check_convergence(ri_l):
                self.convergence_l = True
                break

            self.k_values_liquid = self._update_k_values(self.k_values_liquid, ri_l)

            if self._check_trivial_solution(self.k_values_liquid):
                self.convergence_trivial_solution_l = True
                break

            i += 1

            if i > 100000:
                raise StopIterationError('Число итераций теста стабильности превысило 100000')

    @property
    def k_vals_for_sat_point_calculation(self):
        if (self.stable is None) or self.stable:
            return None
        else:
            if (self.S_v > 1) and (self.S_l > 1):
                if self.S_v > self.S_l:
                    k_vals = self.k_values_vapour
                else:
                    k_vals = self.k_values_liquid
            else:
                k_vals = self.k_values_for_flash

            return k_vals
