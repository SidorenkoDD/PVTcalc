"""
Центральная модель состава флюида.

`Composition` — не dataclass, а "колоночная" структура: мольные доли
(`composition: Dict[str, float]`) + предвычисленные свойства по каждому
компоненту (`composition_data: Dict[property_name, Dict[component_name, value]]`).
Статические свойства чистых компонентов читаются один раз при импорте модуля
из `Utils/DB.json` (глобальная переменная `db`); свойства псевдокомпонентов
C7+ считаются на лету корреляциями из `PlusComponents/`.

Типичный порядок использования::

    comp = Composition({'C1': 0.7, 'C2': 0.2, 'C3': 0.1}, T_res=350.0)
    comp.evaluate_composition_data(c7_plus_correlations={...})  # Tc/Pc/ω/Vc/shift/Kw + параметры EOS
    # comp теперь готов для EOS.BrusilovskiyEOS / PhaseStability / VLE

Это входная точка для всего расчётного пайплайна (см. CLAUDE.md, раздел
Calculation Flow) и то, что дальше передаётся в `Flash`, `TwoPhaseStabilityTest`,
`BrusilovskiyEOS` и эксперименты (`DLE`/`CCE`/`SeparatorTest`).
"""

import logging
from copy import deepcopy

from calc_core.EOS.BaseEOS import EOSType
from calc_core.PlusComponents.PlusComponentCorrelations import PlusComponentProperties
from calc_core.Utils import BRS_EOS_DB as BRSDB
from calc_core.Utils.Errors import InvalidMolarFractionError, NoComponentError
from calc_core.Utils.JsonDBReader import JsonDBReader

logger = logging.getLogger(__name__)
db = JsonDBReader().load_database('DB.json')


def check_and_sort_composition(composition_dict: dict):
    """
    Валидирует состав против БД компонентов и сортирует его в канонический
    порядок: сначала некарбоновые компоненты (N2, CO2, H2S и т.п.) по
    `sequence_number`, затем углеводородные — тоже по `sequence_number`.

    Parameters
    ----------
    composition_dict : dict
        `{имя_компонента: мольная_доля}`. Имена должны присутствовать
        в `db['available_components']` (`Utils/DB.json`).

    Returns
    -------
    dict
        Тот же состав, отсортированный по канону БД.

    Raises
    ------
    NoComponentError
        Если компонент отсутствует в БД.
    """
    # db = JsonDBReader().load_database('DB.json')
    available_components = db['available_components']
    sequence_numbers = db['sequence_number']
    carbon_flags = db['carbon_flag']

    normalized_composition = {}
    c_composition = {}
    nc_composition = {}
    for component in composition_dict:

        if component not in available_components:
            raise NoComponentError(f'Компонент {component} отсутсвует в БД')

        molar_fraction = composition_dict[component]

        seqnum = sequence_numbers[component]
        carbonflag = bool(carbon_flags[component])

        if carbonflag:
            c_composition[(component, seqnum)] = molar_fraction
        else:
            nc_composition[(component, seqnum)] = molar_fraction

    c_composition = dict(sorted(c_composition.items(), key=lambda item: item[0][1]))
    nc_composition = dict(sorted(nc_composition.items(), key=lambda item: item[0][1]))

    for key in nc_composition:
        normalized_composition[key[0]] = nc_composition[key]
    for key in c_composition:
        normalized_composition[key[0]] = c_composition[key]

    return normalized_composition

def _normalize_composition(composition_dict: dict):
    """
    Нормализует мольные доли так, чтобы их сумма была равна 1.

    Parameters
    ----------
    composition_dict : dict
        `{имя_компонента: мольная_доля}` (доли не обязаны быть уже
        нормированы — функция делит каждую на сумму всех).

    Returns
    -------
    dict
        Тот же набор компонент с нормированными долями.

    Raises
    ------
    InvalidMolarFractionError
        Если какая-то мольная доля <= 0.
    """
    normalized_composition = {}

    mol_sum = sum(list(composition_dict.values()))
    for component in composition_dict:
        molar_fraction = composition_dict[component]

        if molar_fraction <= 0.0:
            raise InvalidMolarFractionError(
                f'Мольная доля компонента {component} д.б. больше 0. Передано: {molar_fraction}')

        normalized_composition[component] = molar_fraction / mol_sum
    return normalized_composition


class Composition:
    """
    Состав флюида + свойства компонент, необходимые для EOS.

    Хранит мольные доли (`.composition`) и "колоночную" структуру свойств
    (`.composition_data`) — по одному словарю `{компонент: значение}` на
    каждое свойство (молярная масса, Tc, Pc, ацентрический фактор, параметры
    a/b/c/d EOS Брусиловского, матрица BIP и т.д.). Объекта отдельного
    компонента в живом коде нет — см. CLAUDE.md.
    """

    def __init__(self, zi: dict, / , T_res: float = 373.15, eos_name: EOSType = EOSType.PREOS):
        """
        Parameters
        ----------
        zi : dict
            `{имя_компонента: мольная_доля}`. Валидируется и сортируется
            через `check_and_sort_composition` (компоненты должны быть
            в `Utils/DB.json`). Доли не обязаны суммироваться в 1 на этом
            шаге — см. `normalize_composition()`.
        T_res : float, optional
            Пластовая (расчётная) температура, K. По умолчанию 373.15.
        eos_name : EOSType, optional
            Тип EOS — влияет на выбор BIP-таблицы и shift-параметра
            (`PREOS`/`SRKEOS` берут готовые значения из БД, `BRSEOS`
            считает их корреляциями в `_calc_bip_brs`). По умолчанию `PREOS`.

        После инициализации `composition_data` заполнен только частично
        (молярная масса, Tb, gamma, флаги C5+/C7+/carbon — из
        `_fill_M_Tb_gamma_from_db`); полный расчёт свойств делает
        `evaluate_composition_data(...)`.
        """

        self._composition = check_and_sort_composition(zi)
        self._T_res = T_res  # Пластовая температура
        self._eos_name = eos_name  # Уравнение состояния

        self._composition_data = {
            'molar_mass': {},
            'gamma': {},
            'Tb': {},
            'critical_pressure': {},
            'critical_temperature': {},
            'acentric_factor': {},
            'critical_volume': {},
            'shift_parameter': {},
            'peneloux_correction': {},
            'bip': {},
            'dbipdT': {},
            'c5_plus_flag': {},
            'c7_plus_flag': {},
            'carbon_flag': {},
            'Kw': {},
            'Zc': {},
            'OmegaC': {},
            'Psi': {},
            'alpha': {},
            'beta': {},
            'sigma': {},
            'delta': {},
            'a': {},
            'b': {},
            'c': {},
            'd': {},
        }

        self._fill_M_Tb_gamma_from_db()

    def _evaluate_bips(self, only_brs=False):
        """
        Заполняет `composition_data['bip']`/`['dbipdT']` — полную (nc x nc)
        матрицу коэффициентов парного взаимодействия для текущего `eos_name`.

        Parameters
        ----------
        only_brs : bool, optional
            Если True — пересчитать BIP только если `eos_name == BRSEOS`
            (BIP для Брусиловского зависят от температуры, поэтому
            пересчитываются при каждой смене `.T`; BIP для PR/SRK от
            температуры не зависят и повторно не считаются). По умолчанию False.
        """
        if not only_brs:
            if self._eos_name == EOSType.PREOS:
                bip_func = self._calc_bip_pr
            elif self._eos_name == EOSType.SRKEOS:
                bip_func = self._calc_bip_srk
            else: # self._eos_name == 'BRSEOS'
                bip_func = self._calc_bip_brs

            for comp1 in list(self._composition):
                for comp2 in list(self._composition):
                    # self._composition_data['bip'][comp1][comp2] = bip_func(comp1, comp2, db)
                    bip, dbipdT = bip_func(comp1, comp2)
                    self._composition_data['bip'].setdefault(comp1, {})[comp2] = bip
                    self._composition_data['dbipdT'].setdefault(comp1, {})[comp2] = dbipdT

        elif self._eos_name == EOSType.BRSEOS:
            bip_func = self._calc_bip_brs
            for comp1 in list(self._composition):
                for comp2 in list(self._composition):
                    # self._composition_data['bip'][comp1][comp2] = bip_func(comp1, comp2, db)
                    bip, dbipdT = bip_func(comp1, comp2)
                    self._composition_data['bip'].setdefault(comp1, {})[comp2] = bip
                    self._composition_data['dbipdT'].setdefault(comp1, {})[comp2] = dbipdT

        else:
            pass

    def _calc_bip_pr(self, component1: str, component2: str):
        """
        BIP для Пенга-Робинсона: из БД (`db['bip_pr']`), либо 0.0, если
        хотя бы один из компонентов — C7+ (TODO в коде: для C7+ пар
        значение не считается, а просто зануляется).

        Parameters
        ----------
        component1, component2 : str

        Returns
        -------
        tuple[float, float]
            `(bip, dbipdT)` — второе значение всегда 0.0 (BIP для PR
            от температуры не зависит).
        """
        c7_plus_flag1 = self._composition_data['c7_plus_flag'][component1]
        c7_plus_flag2 = self._composition_data['c7_plus_flag'][component2]

        if any([c7_plus_flag1, c7_plus_flag2]):
            # TODO
            bip = 0.0
        else:
            bip = db['bip_pr'][component1][component2]
        return bip, 0.0

    def _calc_bip_srk(self, component1: str, component2: str):
        """
        BIP для Соаве-Редлиха-Квонга: из БД (`db['bip_srk']`), либо 0.0, если
        хотя бы один из компонентов — C7+.

        Parameters
        ----------
        component1, component2 : str

        Returns
        -------
        tuple[float, float]
            `(bip, dbipdT)` — второе значение всегда 0.0.
        """
        c7_plus_flag1 = self._composition_data['c7_plus_flag'][component1]
        c7_plus_flag2 = self._composition_data['c7_plus_flag'][component2]

        if any([c7_plus_flag1, c7_plus_flag2]):
            bip = 0.0
        else:
            bip = db['bip_srk'][component1][component2]
        return bip, 0.0

    def _calc_bip_brs(self, component1: str, component2: str):
        """
        BIP для EOS Брусиловского. Если хотя бы один компонент из фракции C5+ —
        считает коэффициент корреляцией (`BRS_EOS_DB.evaluate_BRS_EOS_bip_for_c5_plus`,
        зависит от температуры, поэтому и от `dbipdT`); иначе — по таблице
        коэффициентов для лёгких компонент (`evaluate_BRS_EOS_bip_below_c5_plus`).

        Parameters
        ----------
        component1, component2 : str

        Returns
        -------
        tuple[float, float]
            `(bip, dbipdT)`.
        """
        c5_plus_flag1 = self._composition_data['c5_plus_flag'][component1]
        c5_plus_flag2 = self._composition_data['c5_plus_flag'][component2]

        if any([c5_plus_flag1, c5_plus_flag2]):
            # comp1 - Из фракции C5+
            (comp1, _), (comp2, _) = list(sorted([(component1, c5_plus_flag1), (component2, c5_plus_flag2)],
                                                 key=lambda x: x[1], reverse=True))

            omega = self._composition_data['acentric_factor'][comp1]
            Kw = self._composition_data['Kw'][comp1]
            # logger.debug(f'{comp1} {omega} {Kw}')

            return BRSDB.evaluate_BRS_EOS_bip_for_c5_plus(comp1, comp2, self._T_res, omega, Kw)
        else:
            return BRSDB.evaluate_BRS_EOS_bip_below_c5_plus(component1, component2, self._T_res)

    def _fill_M_Tb_gamma_from_db(self):
        """
        Заполняет базовые свойства (`molar_mass`, `gamma`, `Tb`, флаги
        `c7_plus_flag`/`c5_plus_flag`/`carbon_flag`) из `Utils/DB.json`
        для всех компонент состава. Вызывается автоматически в `__init__`.
        """
        for component in self._composition:
            self._composition_data['molar_mass'][component] = db['molar_mass'][component]
            self._composition_data['gamma'][component] = db['gamma'][component]
            self._composition_data['Tb'][component] = db['Tb'][component]
            self._composition_data['c7_plus_flag'][component] = bool(db['c7_plus_flag'][component])
            self._composition_data['c5_plus_flag'][component] = bool(db['c5_plus_flag'][component])
            self._composition_data['carbon_flag'][component] = bool(db['carbon_flag'][component])

    def evaluate_composition_data(self, c7_plus_correlations: dict):
        """
        Главная точка входа: полный расчёт `composition_data` для всех
        компонент состава. Вызывать один раз после `__init__`, прежде чем
        передавать `Composition` в `BrusilovskiyEOS`/`Flash`/etc.

        Для C7+ компонент (`c7_plus_flag`) свойства (Tc, Pc, ω, Vc, shift,
        Kw) считаются корреляциями через `PlusComponentProperties`
        (см. `c7_plus_correlations`); для обычных компонент — берутся из БД
        напрямую (с учётом `eos_name` для shift-параметра: PR/SRK берут
        готовое значение из БД, BRS получает 0.0 здесь и настоящий shift
        считается позже через `peneloux_correction`).

        После этого для каждой компоненты вызывается
        `_evaluate_eos_dependent_props_for_component` (параметры a/b/c/d
        EOS) и заполняется вся матрица BIP (`_evaluate_bips`).

        Parameters
        ----------
        c7_plus_correlations : dict
            Конфиг выбора корреляций для C7+, например
            `{'critical_temperature': 'pedersen', 'critical_pressure': 'rizari_daubert',
            'acentric_factor': 'riazi al sahhaf', 'critical_volume': 'hall yarborough',
            'Kw': 'k watson', 'shift_parameter': 'jhaveri youngren'}`
            (см. `PlusComponents/PlusComponentCorrelations.py` за списком
            допустимых ключей корреляций). Может быть пустым словарём только
            если в составе вообще нет C7+ компонент — иначе упадёт `ValueError`
            внутри `PlusComponentProperties`.
        """

        # Расчет свойств отдельных компонент
        for component in self._composition:
            c7_plus_flag = self._composition_data['c7_plus_flag'][component]
            c5_plus_flag = self._composition_data['c5_plus_flag'][component]
            M = self._composition_data['molar_mass'][component]
            gamma = self._composition_data['gamma'][component]
            Tb = self._composition_data['Tb'][component]

            if c7_plus_flag:
                prop_calculator = PlusComponentProperties(gamma=gamma, Tb=Tb, M=M,
                                                          correlations_config=c7_plus_correlations)
                prop_calculator.calculate_all()

                Kw = prop_calculator.data['Kw']
                Tc = prop_calculator.data['critical_temperature']
                Pc = prop_calculator.data['critical_pressure']
                af = prop_calculator.data['acentric_factor']
                Vc = prop_calculator.data['critical_volume']
                Sshift = prop_calculator.data['shift_parameter']

            else:
                if c5_plus_flag:
                    Kw = db['Kw'][component]
                else:
                    Kw = None
                Pc = db['critical_pressure'][component]
                Tc = db['critical_temperature'][component]
                af = db['acentric_factor'][component]
                Vc = db['critical_volume'][component]

                if self._eos_name == EOSType.PREOS:
                    Sshift = db['shift_parameter_pr'][component]
                elif self._eos_name == EOSType.SRKEOS:
                    Sshift = db['shift_parameter_srk'][component]
                else:
                    Sshift = 0.0

            self._composition_data['Kw'][component] = Kw
            self._composition_data['critical_pressure'][component] = Pc
            self._composition_data['critical_temperature'][component] = Tc
            self._composition_data['acentric_factor'][component] = af
            self._composition_data['critical_volume'][component] = Vc
            self._composition_data['shift_parameter'][component] = Sshift

            # Расчет свойств, зависящих от УРС и температуры
            self._evaluate_eos_dependent_props_for_component(component)

        # Расчет коэффициентов парного взаимодействия
        #  TODO: Если значения нет в БД - расчет по Tb (см. Брусиловский)) (Расширить БД bip_pr до C10)
        self._evaluate_bips()

    def _evaluate_eos_dependent_props_for_component(self, component):
        """
        Считает параметры EOS Брусиловского (OmegaC/Zc/Psi → alpha/beta/sigma/delta
        → a/b/c/d, через `Utils/BRS_EOS_DB.py`) и `peneloux_correction`
        (= b * shift_parameter) для одной компоненты. Требует, чтобы для неё
        уже были заполнены `acentric_factor`/`critical_pressure`/`critical_temperature`
        (иначе просто логирует ошибку и ничего не считает — не бросает исключение).

        Вызывается из `evaluate_composition_data` (для всех компонент) и из
        `edit_component_properties`/`T.setter` (для пересчёта после точечного
        изменения свойства или температуры).

        Parameters
        ----------
        component : str
        """

        af = self._composition_data['acentric_factor'].get(component)
        c5_plus_flag = self._composition_data['c5_plus_flag'].get(component)
        Pc = self._composition_data['critical_pressure'].get(component)
        Tc = self._composition_data['critical_temperature'].get(component)
        sshift = self._composition_data['shift_parameter'].get(component, 0.0)

        if af is None:
            logger.error(f'Отсутствует рассчитанное значение ацентрического фактора для компонента: {component}')
        elif Pc is None:
            logger.error(f'Отсутствует рассчитанное значение критического давления для компонента: {component}')
        elif Tc is None:
            logger.error(f'Отсутствует рассчитанное значение критической температуры для компонента: {component}')
        else:
            OmegaC, Zc, Psi = BRSDB.evaluate_BRS_EOS_component_params(component, af, self._eos_name, c5_plus_flag)
            self._composition_data['OmegaC'][component] = OmegaC
            self._composition_data['Zc'][component] = Zc
            self._composition_data['Psi'][component] = Psi

            alpha, beta, sigma, delta = BRSDB.calc_alpha_beta_gamma_delta(OmegaC, Zc)
            self._composition_data['alpha'][component] = alpha
            self._composition_data['beta'][component] = beta
            self._composition_data['sigma'][component] = sigma
            self._composition_data['delta'][component] = delta

            a, b, c, d = BRSDB.calc_a_b_c_d(self._T_res, Pc, Tc, alpha, beta, sigma, delta, Psi)
            self._composition_data['a'][component] = a
            self._composition_data['b'][component] = b
            self._composition_data['c'][component] = c
            self._composition_data['d'][component] = d

            self._composition_data['peneloux_correction'][component] = b * sshift

    def _evaluate_eos_dependent_props(self):
        """
        Пересчитывает EOS-зависимые свойства (см.
        `_evaluate_eos_dependent_props_for_component`) для *всех* компонент
        состава разом. Вызывается при смене `.T` или `.eos_name`.
        """
        for component in self._composition:
            self._evaluate_eos_dependent_props_for_component(component)

    def recalculate_property(self):
        """Заготовка на будущее — пересчёт отдельного свойства по корреляции. Не реализовано (TODO в коде)."""
        # TODO: Для перерасчета отдельного свойства по корреляции
        pass

    def new_composition(self, molar_fractions: dict, deep_copy=False):
        """
        Возвращает копию состава с переданными мольными долями (компонентный состав должен совпадать!)

        Основной механизм передачи состава жидкости/пара между ступенями в
        `TwoPhaseStabilityTest`/`PhaseEquilibriumNewton` (для каждой пробной
        фазы) и между стадиями PVT-экспериментов (`DLE`/`CCE`/`SeparatorTest`).
        Новый набор долей нормализуется (`_normalize_composition`); свойства
        компонент (`composition_data`) переиспользуются из текущего объекта —
        по ссылке (`deep_copy=False`, дефолт) или полной копией (`deep_copy=True`).

        Parameters
        ----------
        molar_fractions : dict
            `{имя_компонента: мольная_доля}` — набор компонент должен
            **точно совпадать** с текущим составом (не обязательно по
            значениям, но по множеству имён).
        deep_copy : bool, optional
            `False` (по умолчанию) — новый объект разделяет `composition_data`
            с исходным (экономит память, но правки одного видны в другом).
            `True` — независимая глубокая копия.

        Returns
        -------
        Composition
            Новый объект с теми же `T_res`/`eos_name`, что у исходного.

        Raises
        ------
        RuntimeError
            Если набор компонент не совпадает с исходным.
        """
        new_fractions = check_and_sort_composition(molar_fractions)
        if not ((len(new_fractions) == len(self._composition))
                and (set(list(new_fractions.keys())) == set(self._composition.keys()))):
            raise RuntimeError('Переданный компонентный состав не идентичен исходному')
        new_fractions = _normalize_composition(new_fractions)
        new_composition = Composition(new_fractions, T_res=self._T_res, eos_name=self._eos_name)
        new_composition._composition_data = self._composition_data if deep_copy is False else deepcopy(self._composition_data)
        return new_composition

    def normalize_composition(self):
        """Нормализует `self.composition` так, чтобы мольные доли суммировались в 1 (см. `_normalize_composition`)."""
        self._composition = _normalize_composition(self._composition)

    def edit_component_properties(self, component: str, property_name: str, value: float):
        """
        Метод позволяет изменять значение свойства для компонента

        :param component: Название компонента
        :param property_name: Название свойства ('molar_mass', 'gamma', 'Tb', 'critical_pressure',
        'critical_temperature', 'acentric_factor', 'shift_parameter', 'critical_volume')
        :param value: Значение свойства
        """

        avalible_properties = ['molar_mass', 'gamma', 'Tb', 'critical_pressure', 'critical_temperature',
                               'acentric_factor', 'shift_parameter', 'critical_volume', 'Kw', 'peneloux_correction']

        if property_name not in avalible_properties:
            raise KeyError(f'Свойство {property_name} недоступно для изменения или отсутствует!')

        if property_name == 'peneloux_correction':
            self.edit_cpen_for_component(component, value)
        else:
            self._composition_data[property_name][component] = value

        if property_name in ['critical_pressure', 'critical_temperature', 'acentric_factor']:
            self._evaluate_eos_dependent_props_for_component(component)

    def edit_bip_for_components(self, component1: str, component2: str, bip: float):
        """
        Точечно перезаписывает коэффициент парного взаимодействия для пары
        компонент (симметрично — `bip[c1][c2]` и `bip[c2][c1]`).

        Parameters
        ----------
        component1, component2 : str
            Должны присутствовать в текущем составе.
        bip : float
            Новое значение BIP.

        Raises
        ------
        KeyError
            Если один из компонентов отсутствует в составе.
        """
        if component1 not in self._composition:
            raise KeyError(f'Компонент {component1} отсутствует в составе!')
        if component2 not in self._composition:
            raise KeyError(f'Компонент {component2} отсутствует в составе!')
        self._composition_data['bip'][component1][component2] = bip
        self._composition_data['bip'][component2][component1] = bip

    def edit_dbipdt_for_components(self, component1: str, component2: str, dbipdt: float):
        """
        Точечно перезаписывает производную BIP по температуре для пары
        компонент (симметрично).

        Parameters
        ----------
        component1, component2 : str
            Должны присутствовать в текущем составе.
        dbipdt : float

        Raises
        ------
        KeyError
            Если один из компонентов отсутствует в составе.
        """
        if component1 not in self._composition:
            raise KeyError(f'Компонент {component1} отсутствует в составе!')
        if component2 not in self._composition:
            raise KeyError(f'Компонент {component2} отсутствует в составе!')
        self._composition_data['dbipdT'][component1][component2] = dbipdt
        self._composition_data['dbipdT'][component2][component1] = dbipdt

    def edit_cpen_for_component(self, component: str, cpen: float):
        """
        Задаёт `peneloux_correction` напрямую (в единицах объёма) и обратно
        пересчитывает соответствующий `shift_parameter` (= cpen / b).
        Вызывается из `edit_component_properties(..., 'peneloux_correction', ...)`.

        Parameters
        ----------
        component : str
        cpen : float | None
            Новое значение поправки Пенелу. `None` — сбросить и shift, и
            поправку в 0.0.
        """
        b = self._composition_data['b'].get(component)
        if cpen is not None:
            self._composition_data['shift_parameter'][component] = cpen / b
            self._composition_data['peneloux_correction'][component] = cpen
        else:
            self._composition_data['shift_parameter'][component] = 0.0
            self._composition_data['peneloux_correction'][component] = 0.0

    @property
    def composition(self):
        """dict: `{имя_компонента: мольная_доля}` — возвращается по ссылке, не копией."""
        return self._composition
        # return deepcopy(self._composition)

    @property
    def composition_data(self):
        """dict: `{property_name: {имя_компонента: значение}}` — возвращается по ссылке, не копией."""
        return self._composition_data
        # return deepcopy(self._composition_data)

    @property
    def T(self):
        """float: текущая температура состава, K."""
        return self._T_res

    @T.setter
    def T(self, T: float):
        """
        Смена температуры **пересчитывает** EOS-зависимые свойства всех
        компонент (`_evaluate_eos_dependent_props`) и BIP для Брусиловского
        (`_evaluate_bips(only_brs=True)`) — не побочный эффект, а основной
        механизм актуализации состояния при переходе на новую P/T-точку
        (так вызывается из `Flash.__init__`).
        """
        self._T_res = float(T)
        self._evaluate_eos_dependent_props()
        self._evaluate_bips(only_brs=True)

    @property
    def eos_name(self):
        """EOSType: текущий тип EOS (`PREOS`/`SRKEOS`/`BRSEOS`), влияет на выбор BIP и shift-параметра."""
        return self._eos_name

    @eos_name.setter
    def eos_name(self, eos_name: EOSType):
        """Смена типа EOS пересчитывает EOS-зависимые свойства и всю матрицу BIP заново."""
        self._eos_name = eos_name
        self._evaluate_eos_dependent_props()
        self._evaluate_bips()


    @classmethod
    def from_db(cls, db_path: str = "models.json", T_res: float = 373.15):
        """
        Загружает сохранённые "снэпшоты" моделей флюидов из JSON-файла
        (обычно `models.json` в корне репозитория, пишется через
        `Utils/Export.py::ModelJSONDB`) и возвращает объект-прокси, через
        который каждая модель доступна как атрибут по её ключу в файле.

        Composition_data для каждой модели инжектируется напрямую из JSON
        (`obj._composition_data = rec["composition_data"]`), **без**
        пересчёта корреляций — то есть `evaluate_composition_data(...)`
        вызывать не нужно (и не стоит, если хотите именно сохранённые числа).

        Parameters
        ----------
        db_path : str, optional
            Путь к JSON-файлу со снэпшотами. По умолчанию `"models.json"`.
        T_res : float, optional
            Пластовая температура, которая будет проставлена каждому
            загруженному объекту (в файле не хранится — `Export.py` её
            не сохраняет). По умолчанию 373.15.

        Returns
        -------
        _CompositionDBProxy
            `db.<ключ_модели>` -> `Composition`; `db.list_models()` -> список
            доступных ключей; `dir(db)` тоже их показывает (для автокомплита в IDE).

        Raises
        ------
        FileNotFoundError
            Если `db_path` не существует.

        Пример
        ------
        >>> db = Composition.from_db("models.json")
        >>> comp = db.KRSNL_PVTSIM
        """
        from pathlib import Path

        from calc_core.Utils.ModelStore import read_model_store

        path = Path(db_path)
        if not path.exists():
            raise FileNotFoundError(f"Файл БД не найден: {path}")

        db_data = read_model_store(path)

        class _CompositionDBProxy:
            """Ленивый доступ к моделям из `models.json` через атрибуты (см. `Composition.from_db`)."""

            def __init__(self, data, comp_cls, default_T):
                self._data = data
                self._comp_cls = comp_cls
                self._default_T = default_T

            def __getattr__(self, name: str):
                """Строит `Composition` из записи `self._data[name]` при первом обращении (не кэширует)."""
                if name.startswith('_'):
                    raise AttributeError(name)
                if name not in self._data:
                    raise AttributeError(f"Состав '{name}' не найден в БД. Доступные: {list(self._data.keys())}")

                rec = self._data[name]
                zi = rec["composition"]

                # 1. Стандартная инициализация (создает пустую структуру _composition_data и валидирует zi)
                obj = self._comp_cls(zi)

                # 2. Мгновенный инжект сохраненных свойств из JSON (без пересчета корреляций)
                if "composition_data" in rec:
                    obj._composition_data = rec["composition_data"]

                # 3. Восстановление типа УРС из строки/Enum
                eos_raw = rec.get("eos", "PREOS")
                eos_str = str(eos_raw).upper()
                if "SRK" in eos_str:
                    obj._eos_name = EOSType.SRKEOS
                elif "BRS" in eos_str:
                    obj._eos_name = EOSType.BRSEOS
                else:
                    obj._eos_name = EOSType.PREOS

                # 4. Пластовая температура: новые записи могут хранить её в
                #    ключе "T_res" (см. Export.py, параметр t_res); старые
                #    записи ключа не имеют — используем default_T, как раньше.
                obj._T_res = float(rec.get("T_res", self._default_T))

                return obj

            def __dir__(self):
                """dir(proxy) показывает все ключи моделей — подсказка для автокомплита в IDE."""
                # Подсказка для IDE: dir(proxy) покажет все id моделей
                return list(self._data.keys()) + list(super().__dir__())

            def list_models(self) -> list:
                """list[str]: все доступные ключи моделей в загруженном файле."""
                return list(self._data.keys())

        return _CompositionDBProxy(db_data, cls, T_res)
