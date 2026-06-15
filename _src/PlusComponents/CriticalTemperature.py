import math
from typing import Dict, Callable
from _src.Utils.Constants import CONSTANT_R

class CriticalTemperatureCorrelation:
    """
    Класс, содержащий корреляции для расчета критической температуры
    """

    @staticmethod
    def _roess(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Роесса (Roess, 1936)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическая температура в K
        """
        Tb_F = Tb_K * 9 / 5 - 459.67
        Tc_R = 645.83 + 1.6667 * (gamma * (Tb_F + 100)) - 0.727e-3 * math.pow(gamma * (Tb_F + 100), 2)
        return Tc_R * 5 / 9

    @staticmethod
    def _nokey(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Нокея (Nokey, 1959)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическая температура в K
        """
        Tb_R = Tb_K * 9 / 5
        Tc_R = 19.07871 * math.pow(Tb_R, 0.62164) * math.pow(gamma, 0.2985)
        return Tc_R * 5 / 9

    @staticmethod
    def _cavett(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Каветта (Cavett, 1962)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K
        :return: Критическая температура в K
        """
        gamma_api = 141.5 / gamma - 131.5
        Tb_F = Tb_K * 9 / 5 - 459.67
        Tc_R = (768.07121 + 1.7133693 * Tb_F - 0.10834003 * 1e-2 * math.pow(Tb_F, 2) -
                0.89212579 * 1e-2 * Tb_F * gamma_api +
                0.38890584 * 1e-6 * math.pow(Tb_F, 3) +
                0.5309492 * 1e-5 * math.pow(Tb_F, 2) * gamma_api +
                0.327116 * 1e-7 * math.pow(Tb_F, 2) * math.pow(gamma_api, 2))
        return Tc_R * 5 / 9

    @staticmethod
    def _kesler_lee(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Кеслера-Ли (Kesler-Lee, 1976)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K
        :return: Критическая температура в K
        """
        Tb_R = Tb_K * 9 / 5
        Tc_R = (341.7 + 811 * gamma + (0.4244 + 0.1174 * gamma) * Tb_R + ((0.4669 - 3.2623 * gamma) * 1e5) / Tb_R)
        return Tc_R * 5 / 9

    @staticmethod
    def _pedersen(gamma: float, M: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Педерсен (Pedersen, 1988)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса в г/моль (g/mol)
        :return: Критическая температура в K
        """
        return 163.12 * gamma + 86.052 * math.log(M) + 0.43475 * M - (1877.4 / M)

    @staticmethod
    def _standing(gamma: float, M: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Стэндинга (Standing, 1977)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса в г/моль (g/mol)
        :return: Критическая температура в K
        """
        if M <= 71.2:
            raise ValueError(f'Молярная масса {M} д.б. строго больше минимального значения 71.2')
        return 338 + 202 * math.log10(M - 71.2) + (1361 * math.log10(M) - 2111) * math.log10(gamma)

    @staticmethod
    def _sim_daubert(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Сима-Дауберта (Sim-Daubert, 1980)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическая температура в K
        """
        Tb_R = Tb_K * 9 / 5
        Tc_R = math.exp(3.9934718 * math.pow(Tb_R, 0.08615) * math.pow(gamma, 0.04614))
        return Tc_R * 5 / 9

    @staticmethod
    def _riazi_daubert(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Риази-Дауберта (Riazi-Daubert, 1980)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическая температура в K
        """
        Tb_R = Tb_K * 9 / 5
        Tc_R = 24.27871 * math.pow(Tb_R, 0.58848) * math.pow(gamma, 0.3596)
        return Tc_R * 5 / 9

    @staticmethod
    def _mogoulas_tassios(gamma: float, M: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Могоуласа-Тассиоса (Mogoulas, Tassios, 1990)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса в г/моль (g/mol)

        :return: Критическая температура, K
        """

        Tc_R = -1274.4 + 0.792 * M + 1971 * gamma - 27000 / M + 707.4 / gamma
        return Tc_R * 5 / 9

    @staticmethod
    def _twu(Tb_K: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Тву (Twu, 1984)
        :param Tb_K: Температура кипения в K

        :return: Критическая температура в K
        """
        Tb_R = Tb_K * 9 / 5
        Tc_R = Tb_R / (0.533272 + 0.191017 * 1e-3 * Tb_R + 0.779681 * 1e-7 * math.pow(Tb_R, 2)
                       - 0.284376 * 1e-10 * math.pow(Tb_R, 3) + 0.959468 * 1e28 / math.pow(Tb_R, 13))
        return Tc_R * 5 / 9

    @staticmethod
    def _watansiri_owens_starling(gamma: float, Tb_K: float, M: float) -> float:
        """
        Рассчитывает критическую температуру по корреляции Ватансири-Овенса-Старлинга (Watansiri, Owens, Starling, 1985)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K
        :param M: Молярная масса в г/моль (g/mol)

        :return: Критическая температура в K
        """
        Tb_R = Tb_K * 9 / 5
        ln_Tc_R = (-0.0650504 - 0.0005217 * Tb_R + 0.03095 * math.log(M) + 1.11067 * math.log(Tb_R)
                   + M * (0.078154 * math.pow(gamma, 0.5) - 0.061061 * math.pow(gamma, 1 / 3) - 0.016943 * gamma))
        Tc_R = math.exp(ln_Tc_R)
        return Tc_R * 5 / 9

    def get_correlation(self, method: str) -> Callable:
        """
        Возвращает функцию корреляции по названию
        :param method: Название метода (e.g. Roess, Nokey, Cavett, Kesler Lee, Pedersen, Standing, Sim Daubert, Riazi Daubert, Mogoulas Tassios, Twu, Watansiri Owens Starling)
        """

        method_map = {
            'roess': self._roess,
            'nokey': self._nokey,
            'cavett': self._cavett,
            'kesler lee': self._kesler_lee,
            'pedersen': self._pedersen,
            'standing': self._standing,
            'sim daubert': self._sim_daubert,
            'riazi daubert': self._riazi_daubert,
            'mogoulas tassios': self._mogoulas_tassios,
            'twu': self._twu,
            'watansiri owens starling': self._watansiri_owens_starling,
        }
        _method = method.lower()
        if not _method in method_map:
            raise ValueError(f"CriticalTemperatureCorrelation: Unknown correlation method: {method}")
        return method_map[_method]

    @staticmethod
    def get_required_params(method: str) -> list:
        """
        Возвращает список требуемых параметров для корреляции
        """
        params_map = {
            'roess': ['gamma', 'Tb_K'],
            'nokey': ['gamma', 'Tb_K'],
            'cavett': ['gamma', 'Tb_K'],
            'kesler lee': ['gamma', 'Tb_K'],
            'pedersen': ['gamma', 'M'],
            'standing': ['gamma', 'M'],
            'sim daubert': ['gamma', 'Tb_K'],
            'riazi daubert': ['gamma', 'Tb_K'],
            'mogoulas tassios': ['gamma', 'M'],
            'twu': ['Tb_K', ],
            'watansiri owens starling': ['gamma', 'Tb_K', 'M'],
        }
        _method = method.lower()
        if _method not in params_map:
            raise ValueError(f"CriticalTemperatureCorrelation: Unknown correlation method: {method}")
        return params_map[_method]
