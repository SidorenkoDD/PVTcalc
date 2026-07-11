"""Корреляции ацентрического фактора C7+ (4 метода) — см. `PlusComponentCorrelations.py` за диспетчеризацией."""

import math
from typing import Dict, Callable
from calc_core.Utils.Constants import CONSTANT_R


class AcentricFactorCorrelation:
    """
    Класс, содержащий корреляции для расчета ацентрического фактора
    """

    def get_correlation(self, method: str) -> Callable:
        """
        Возвращает функцию корреляции по названию
        :param method: Название метода (e.g. Edmister, Riazi Al Sahhaf, Mogoulas Tassios, Kesler Lee)
        """
        method_map = {
            'kesler lee': self._kesler_lee,
            'riazi al sahhaf': self._riazi_al_sahhaf,
            'edmister': self._edmister,
            'mogoulas tassios': self._mogoulas_tassios,
        }

        _method = method.lower()
        if not _method in method_map:
            raise ValueError(f"AcentricFactorCorrelation: Unknown correlation method: {method}")
        return method_map[_method]

    @staticmethod
    def get_required_params(method: str) -> list:
        """
        Возвращает список требуемых параметров для корреляции
        :param method: Название метода (e.g. Edmister, Riazi Al Sahhaf, Mogoulas Tassios, Kesler Lee)
        """
        params_map = {
            'edmister': ['Pc_bar', 'Tc_K', 'Tb_K'],
            'riazi al sahhaf': ['M'],
            'kesler lee': ['Pc_bar', 'Tc_K', 'Tb_K', 'Kw'],
            'mogoulas tassios': ['gamma', 'M']
        }
        _method = method.lower()
        if _method not in params_map:
            raise ValueError(f"AcentricFactorCorrelation: Unknown correlation method: {method}")
        return params_map[_method]

    @staticmethod
    def _edmister(Pc_bar: float, Tc_K: float, Tb_K: float) -> float:
        """
        Рассчитывает ацентрический фактор по корреляции Эдмистера (Edmister, 1958)
        :param Pc_bar: Критическое давление, bar
        :param Tc_K: Критическая температура, K
        :param Tb_K: Температура кипения, K

        :return: Ацентрический фактор, []
        """
        Pc_psia = Pc_bar * 14.50377
        Tc_R = Tc_K * 1.8
        Tb_R = Tb_K * 1.8

        return (3 / 7) * (math.log10(Pc_psia / 14.7) / (Tc_R / Tb_R - 1)) - 1

    @staticmethod
    def _riazi_al_sahhaf(M: float) -> float:
        """
        Рассчитывает ацентрический фактор по корреляции Риази - Аль-Саххаф (Riazi, Al-Sahhaf, 1996)
        :param M: Молярная масса в г/моль (g/mol)

        :return: Ацентрический фактор, []
        """
        return - (0.3 - math.exp(-6.252 + 3.64457 * math.pow(M, 0.1)))

    @staticmethod
    def _mogoulas_tassios(gamma: float, M: float) -> float:
        """
        Рассчитывает ацентрический фактор по корреляции Могоуласа-Тассиоса (Mogoulas, Tassios, 1990)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса в г/моль (g/mol)

        :return: Ацентрический фактор, []
        """
        return -0.64235 + 0.00014667 * M + 0.021876 * gamma - 4.559 / M + 0.21699 * math.log(M)

    @staticmethod
    def _kesler_lee(Pc_bar: float, Tc_K: float, Tb_K: float, Kw: float) -> float:
        """
        Рассчитывает ацентрический фактор по корреляции Кеслера-Ли (Kesler, Lee, 1975-1976)
        :param Pc_bar: Критическое давление, bar
        :param Tc_K: Критическая температура, K
        :param Tb_K: Температура кипения, K
        :param Kw: Фактор Ватсона

        :return: Ацентрический фактор, []
        """

        Pc_psia = Pc_bar * 14.50377
        Tc_R = Tc_K * 1.8
        Tb_R = Tb_K * 1.8

        Tbr = Tb_R / Tc_R

        if Tbr < 0.8:
            num = (
                    -math.log(Pc_psia / 14.7)
                    - 5.92714
                    + 6.09648 / Tbr
                    + 1.28862 * math.log(Tbr)
                    - 0.169347 * math.pow(Tbr, 6)
            )
            den = (
                    15.2518
                    - 15.6875 / Tbr
                    - 13.4721 * math.log(Tbr)
                    + 0.43577 * math.pow(Tbr, 6)
            )
            omega = num / den
        else:
            omega = (
                    -7.904
                    + 0.1352 * Kw
                    - 0.007465 * math.pow(Kw, 2)
                    + 8.359 * Tbr
                    + (1.408 - 0.01063 * Kw) / Tbr
            )

        return omega
