import math
from typing import Dict, Callable
from _src.Utils.Constants import CONSTANT_R

class ShiftParameterCorrelation:
    """
    Класс, содержащий корреляции для расчета шифт-параметра
    """

    def get_correlation(self, method: str) -> Callable:
        """
        Возвращает функцию корреляции по названию
        :param method: Название метода (e.g. Jhaveri Youngren, SRK, PR)
        """
        method_map = {
            'jhaveri youngren': self._jhaveri_youngren,
            'srk': self._srk,
            'pr': self._pr,
            '': lambda: 0.0,
            # 'BRS': self._brs,  # TODO
        }

        _method = method.lower()
        if not _method in method_map:
            raise ValueError(f"ShiftParameterCorrelation: Unknown correlation method: {method}")
        return method_map[_method]

    @staticmethod
    def get_required_params(method: str) -> list:
        """
        Возвращает список требуемых параметров для корреляции
        :param method: Название метода (e.g. Jhaveri Youngren, SRK, PR)
        """
        params_map = {
            'jhaveri youngren': ['M', 'Kw'],
            'srk': ['af', ],
            'pr': ['af', ],
            '': [],
            # 'BRS': [],  # TODO

        }
        _method = method.lower()
        if _method not in params_map:
            raise ValueError(f"ShiftParameterCorrelation: Unknown correlation method: {method}")
        return params_map[_method]

    @staticmethod
    def _jhaveri_youngren(M: float, Kw: float) -> float:
        """
        Рассчитывает шифт-параметр s_i (s_i = c_i / b_i) по корреляции Явери-Юнгрен (Jhaveri, Youngren, 1988)
        :param M: Молярная масса, г/моль (g/mol)
        :param Kw: Фактор Ватсона

        :return: Шифт-параметр
        """

        # aromatic
        if 8.5 < Kw < 11:
            a0 = 2.516
            a1 = 0.2008

        # naften
        elif 11 < Kw < 12.5:
            a0 = 3.004
            a1 = 0.2324

        # parafin
        elif 12.5 < Kw < 13.5:
            a0 = 2.258
            a1 = 0.1823

        else:
            return 0.0

        return 1 - (a0 / (math.pow(M, a1)))

    @staticmethod
    def _srk(af: float) -> float:
        """
        Рассчитывает шифт-параметр s_i (s_i = c_i / b_i) для уравнения SRK
        :param af: Ацентрический фактор

        :return: Шифт-параметр
        """
        z_ra = 0.29056 - 0.08775 * af
        # return 0.40768 * (CONSTANT_R * Tc_K / Pc_bar) * (0.29441 - z_ra)
        return 4.70545 * (0.29441 - z_ra)

    @staticmethod
    def _pr(af: float) -> float:
        """
        Рассчитывает шифт-параметр s_i (s_i = c_i / b_i) для уравнения PR
        :param af: Ацентрический фактор

        :return: Шифт-параметр
        """
        z_ra = 0.29056 - 0.08775 * af
        # return 0.50033 * (CONSTANT_R * Tc_K / Pc_bar) * (0.25969 - z_ra)
        return 6.43098 * (0.25969 - z_ra)

    @staticmethod
    def _brs():
        # TODO
        ...
