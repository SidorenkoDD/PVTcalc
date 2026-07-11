import math
from typing import Dict, Callable
from calc_core.Utils.Constants import CONSTANT_R


class CriticalPressureCorrelation:
    """
    Класс, содержащий корреляции для расчета критического давления
    """

    @staticmethod
    def _kesler_lee(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическое давление по корреляции Кеслера-Ли (Kesler, Lee, 1976)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическое давление, bar
        """
        Tb_R = Tb_K * 9 / 5

        first = 8.3634 - 0.0566 / gamma
        second = - (0.24244 + (2.2898 / gamma) + (0.11857 / math.pow(gamma, 2))) * 1e-3 * Tb_R
        third = (1.4685 + (3.648 / gamma) + (0.47227 / math.pow(gamma, 2))) * 1e-7 * math.pow(Tb_R, 2)
        fourth = - (0.42019 + (1.6977 / math.pow(gamma, 2))) * 1e-10 * math.pow(Tb_R, 3)
        ln_Pc_psia = sum([first, second, third, fourth])
        Pc_psia = math.exp(ln_Pc_psia)
        return Pc_psia * 0.0689476

    @staticmethod
    def _riazi_daubert(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическое давление по корреляции Риази-Дауберта (Riazi, Daubert, 1980)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическое давление, bar
        """
        Tb_R = Tb_K * 9 / 5
        Pc_psia = 3.12281 * 1e9 * math.pow(Tb_R, -2.3125) * math.pow(gamma, 2.3201)
        return Pc_psia * 0.0689476

    @staticmethod
    def _cavett(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическое давление по корреляции Каветта (Cavett, 1962)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическое давление, bar
        """
        Tb_F = Tb_K * 9 / 5 - 459.67
        gamma_api = 141.5 / gamma - 131.5

        log_Pc_psia = (
                2.8290406
                + (0.94120109 * 1e-3) * Tb_F
                - (0.30474749 * 1e-5) * math.pow(Tb_F, 2)
                - (0.2087611 * 1e-4) * gamma_api * Tb_F
                + (0.15184103 * 1e-8) * math.pow(Tb_F, 3)
                + (0.11047899 * 1e-7) * gamma_api * math.pow(Tb_F, 2)
                - (0.48271599 * 1e-7) * math.pow(gamma_api, 2) * Tb_F
                + (0.13949619 * 1e-9) * math.pow(gamma_api, 2) * math.pow(Tb_F, 2)
        )

        Pc_psia = math.pow(10, log_Pc_psia)
        return Pc_psia * 0.0689476

    @staticmethod
    def _pedersen(gamma: float, M: float) -> float:
        """
        Рассчитывает критическое давление по корреляции Педерсен (Pedersen, 1988)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса в г/моль (g/mol)

        :return: Критическое давление, bar
        """
        return math.exp(-0.13408 + 2.5019 * gamma + (208.46 / M) - (3987.2 / math.pow(M, 2)))

    @staticmethod
    def _standing(gamma: float, M: float) -> float:
        """
        Рассчитывает критическое давление по корреляции Стэндинга (Standing, 1977)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса в г/моль (g/mol)

        :return: Критическое давление, bar
        """
        if M <= 53.7:
            raise ValueError(f'Молярная масса {M} д.б. больше минимального значения 61.1')
        Pc_mpa = 8.191 - 2.97 * math.log10(M - 61.1) + (15.99 - 5.87 * math.log10(M - 53.7)) * (gamma - 0.8)
        return Pc_mpa * 10

    @staticmethod
    def _sim_daubert(gamma: float, Tb_K: float) -> float:
        """
        Рассчитывает критическое давление по корреляции Сима-Дауберта (Sim, Daubert, 1980)
        :param gamma: Относительная плотность (Specific gravity)
        :param Tb_K: Температура кипения в K

        :return: Критическое давление, bar
        """
        Tb_R = Tb_K * 9 / 5
        Pc_psia = 3.48242 * 1e9 * math.pow(Tb_R, -2.3177) * math.pow(gamma, 2.4853)
        return Pc_psia * 0.0689476

    @staticmethod
    def _mogoulas_tassios(gamma: float, M: float) -> float:
        """
        Рассчитывает критическое давление по корреляции Могоуласа-Тассиоса (Mogoulas, Tassios, 1990)
        :param gamma: Относительная плотность (Specific gravity)
        :param M: Молярная масса в г/моль (g/mol)

        :return: Критическое давление, bar
        """
        ln_Pc_psia = 0.01901 - 0.0048442 * M + 0.13239 * gamma + (227 / M) - (1.1663 / gamma) + 1.2702 * math.log(M)
        Pc_psia = math.exp(ln_Pc_psia)
        return Pc_psia * 0.0689476

    @staticmethod
    def _twu(gamma: float, Tb_K: float):
        # TODO
        ...

    @staticmethod
    def _watansiri_owens_starling(gamma: float, Tb_K: float, M: float) -> float:
        # TODO
        ...

    @staticmethod
    def _pc_from_eos(gamma, M, Tc):
        # TODO
        ...

    def get_correlation(self, method: str) -> Callable:
        """
        Возвращает функцию корреляции по названию
        :param method: Название метода (e.g. Cavett, Kesler Lee, Pedersen, Standing, Sim Daubert, Riazi Daubert, Mogoulas Tassios)
        """
        method_map = {
            'kesler lee': self._kesler_lee,
            'riazi daubert': self._riazi_daubert,
            'cavett': self._cavett,
            'pedersen': self._pedersen,
            'standing': self._standing,
            'sim daubert': self._sim_daubert,
            'mogoulas tassios': self._mogoulas_tassios
        }

        _method = method.lower()
        if not _method in method_map:
            raise ValueError(f"CriticalPressureCorrelation: Unknown correlation method: {method}")
        return method_map[_method]

    @staticmethod
    def get_required_params(method: str) -> list:
        """
        Возвращает список требуемых параметров для корреляции
        :param method: Название метода (e.g. Cavett, Kesler Lee, Pedersen, Standing, Sim Daubert, Riazi Daubert, Mogoulas Tassios)
        """

        params_map = {
            'kesler lee': ['gamma', 'Tb_K'],
            'riazi daubert': ['gamma', 'Tb_K'],
            'cavett': ['gamma', 'Tb_K'],
            'pedersen': ['gamma', 'M'],
            'standing': ['gamma', 'M'],
            'sim daubert': ['gamma', 'Tb_K'],
            'mogoulas tassios': ['gamma', 'M']
        }

        _method = method.lower()
        if _method not in params_map:
            raise ValueError(f"CriticalPressureCorrelation: Unknown correlation method: {method}")
        return params_map[_method]
