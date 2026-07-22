"""Корреляции критического объёма C7+ (4 реализованных метода + 1 заготовка) — см. `PlusComponentCorrelations.py` за диспетчеризацией."""

import math
from typing import Callable

from calc_core.Utils.Constants import CONSTANT_R


class CriticalVolumeCorrelation:
    """
    Класс, содержащий корреляции для расчета критического объема
    """

    def get_correlation(self, method: str) -> Callable:
        """
        Возвращает функцию корреляции по названию
        :param method: Название метода (e.g. Hall Yarborough, Riazi Daubert, Reid, Lohrenz)
        """
        method_map = {
            'hall yarborough': self._hall_yarborough,
            'riazi daubert': self._riazi_daubert,
            'reid': self._reid,
            'lohrenz': self._lohrenz,
        }

        _method = method.lower()
        if _method not in method_map:
            raise ValueError(f"CriticalVolumeCorrelation: Unknown correlation method: {method}")
        return method_map[_method]

    @staticmethod
    def get_required_params(method: str) -> list:
        """
        Возвращает список требуемых параметров для корреляции
        :param method: Название метода (e.g. Hall Yarborough, Riazi Daubert, Reid, Lohrenz)
        """
        params_map = {
            'hall yarborough': ['gamma', 'M'],
            'riazi daubert': ['gamma', 'Tb_K'],
            'reid': ['af', 'Pc_bar', 'Tc_K'],
            'lohrenz': ['gamma', 'M']
        }
        _method = method.lower()
        if _method not in params_map:
            raise ValueError(f"CriticalVolumeCorrelation: Unknown correlation method: {method}")
        return params_map[_method]

    @staticmethod
    def _riazi_daubert(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критический объем по корреляции Риази-Дауберта (Riazi, Daubert, 1980)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критический объем, cm3/mol
        """
        Tb_R = Tb_K * 9 / 5
        Vc_ft3_lbmol = (7.0434 * 1e-7) * math.pow(Tb_R, 2.3829) * math.pow(gamma, -1.683)
        return Vc_ft3_lbmol * 62.42796

    @staticmethod
    def _hall_yarborough(gamma: float, M: float) -> float:
        """
        Рассчитывает критический объем по корреляции Холла-Ярборо (Hall, Yarborough, 1971)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса, г/моль (g/mol)

        :return: Критический объем, cm3/mol
        """
        Vc_ft3_lbmol = 0.025 * math.pow(M, 1.15) * math.pow(gamma, -0.7935)
        return Vc_ft3_lbmol * 62.42796

    @staticmethod
    def _reid(af: float, Pc_bar: float, Tc_K: float) -> float:
        """
        Рассчитывает критический объем по корреляции Рейд (Reid et al., 1977)
        :param af: Ацентрический фактор
        :param Pc_bar: Критическое давление, бар
        :param Tc_K: Критическая температура, K

        :return: Критический объем, cm3/mol
        """
        # CONSTANT_R задан в J/(mol*K). При Pc в bar и результате в cm3/mol
        # нужен коэффициент 10: 1 J = 10 bar*cm3.
        return (0.2918 - 0.0928 * af) * (10.0 * CONSTANT_R) * Tc_K / Pc_bar

    @staticmethod
    def _lohrenz(gamma: float, M: float) -> float:
        """
        Рассчитывает критический объем по корреляции Лоренца (Lohrenz et al., 1964)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса, г/моль (g/mol)

        :return: Критический объем, cm3/mol
        """
        Vc_ft3_lbmol = 21.573 + 0.015122 * M - 27.656 * gamma + 0.070615 * M * gamma
        return Vc_ft3_lbmol * 62.42796

    @staticmethod
    def _riedel():
        """Не реализовано (TODO). Не зарегистрирован в `get_correlation`/`get_required_params` — недостижим через публичный API."""
        # TODO
        ...
