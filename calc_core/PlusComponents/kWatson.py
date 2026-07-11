import math
from typing import Dict, Callable
from calc_core.Utils.Constants import CONSTANT_R

class KWatsonCorrelation:
    """
    Класс, содержащий корреляции для расчета фактора Ватсона
    """

    def get_correlation(self, method: str) -> Callable:
        """
        Возвращает функцию корреляции по названию
        :param method: Название метода (e.g. K Watson, Riazi Daubert)
        """
        method_map = {
            'k watson': self._k_watson,
            'riazi daubert': self._riazi_daubert,
        }

        _method = method.lower()
        if not _method in method_map:
            raise ValueError(f"KWatsonCorrelation: Unknown correlation method: {method}")
        return method_map[_method]

    @staticmethod
    def get_required_params(method: str) -> list:
        """
        Возвращает список требуемых параметров для корреляции
        :param method: Название метода (e.g. K Watson, Riazi Daubert)
        """
        params_map = {
            'k watson': ['gamma', 'Tb_K'],
            'riazi daubert': ['gamma', 'M'],
        }
        _method = method.lower()
        if _method not in params_map:
            raise ValueError(f"KWatsonCorrelation: Unknown correlation method: {method}")
        return params_map[_method]

    @staticmethod
    def _k_watson(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает фактор Ватсона (Universal Oil Products (UOP) Characterization Factor)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Фактор Ватсона
        """
        Tb_R = Tb_K * 9 / 5
        return math.pow(Tb_R, 1 / 3) / gamma

    @staticmethod
    def _riazi_daubert(gamma: float, M: float) -> float:
        """
        Рассчитывает фактор Ватсона (Universal Oil Products (UOP) Characterization Factor)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса, г/моль (g/mol)

        :return: Фактор Ватсона
        """
        return 4.5579 * math.pow(M, 0.15178) * math.pow(gamma, -0.84573)
