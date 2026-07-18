"""
Тонкий фасад над стандартными PVT-экспериментами.

Привязан к `CompositionalModel` (доступен как `model.experiments`) — берёт
состав из `model.composition`, вызывает `DLE`/`CCE`/`SeparatorTest` и логирует
каждый вызов в `model.result_store_object` — по тому же паттерну, что
`PhaseEnvelopeFacade` (`calc_core/PhaseEnvelope/PhaseEnvelopeFacade.py`).

Единообразие, которое даёт этот слой: `DLE.calculate()`/`SeparatorTest.calculate()`
возвращают сырой `list[FlashResult]` (таблица `Bo`/`Rs` — отдельно, в
`self._dle_df`), а `CCE.calculate()` возвращает `pd.DataFrame` напрямую —
здесь все три метода единообразно возвращают `pd.DataFrame`.

Как и в `PhaseEnvelopeFacade`: состав передаётся в `DLE`/`CCE`/`SeparatorTest`
**копией** (`Composition.new_composition(..., deep_copy=True)`), а не общей
ссылкой — иначе `model.composition` мутировался бы (`composition.T = ...`
внутри `DLE.calculate()`/`CCE.__init__` — тяжёлый setter, пересчитывает
EOS-зависимые свойства/BIP). Массивы давлений/температур на входе тоже
копируются (`CCE` дописывает найденный `P_sat` в конец своего
`pressure_arr_bar` на месте).

`CCE.calculate_parallel()` не оборачивается — не проходит постобработку
(`v/v_res` и т.п.) и на практике медленнее последовательного варианта (см.
заметку в `test_notebook.ipynb`, раздел "CCE в параллель"). Явного сигнала
несходимости у `DLE`/`CCE`/`SeparatorTest` нет (в отличие от Bubble/Dew/
Critical point в `PhaseEnvelopeFacade`) — `ConvergenceError` здесь не
поднимается, это вне scope (см. `docs/BACKLOG.md`).
"""

import logging

import pandas as pd

from calc_core.Experiments.DLE import DLE
from calc_core.Experiments.CCE import CCE
from calc_core.Experiments.SeparatorTest import SeparatorTest
from calc_core.Utils.Validation import validate_positive_pressure, validate_temperature_kelvin, validate_temperature_celsius

logger = logging.getLogger(__name__)


class ExperimentsFacade:
    """Доступ к стандартным PVT-экспериментам через `model.experiments`."""

    def __init__(self, model: 'CompositionalModel'):
        """
        Parameters
        ----------
        model : CompositionalModel
            Родительская модель — состав берётся как `model.composition`,
            результаты логируются в `model.result_store_object`.
        """
        self._model = model

    def _composition_copy(self):
        """Независимая копия `model.composition` — см. докстринг модуля."""
        composition = self._model.composition
        return composition.new_composition(composition.composition, deep_copy=True)

    def dle(self, pressure_arr_bar: list, reservoir_pressure_bar: float,
            reservoir_temperature_c: float) -> pd.DataFrame:
        """
        Дифференциальная конденсация (см. `calc_core/Experiments/DLE.py`).

        Parameters
        ----------
        pressure_arr_bar : list
            Ступени давления, бар.
        reservoir_pressure_bar : float
            Пластовое давление, бар.
        reservoir_temperature_c : float
            Пластовая температура, °C.

        Returns
        -------
        pandas.DataFrame
            Таблица по ступеням давления с колонками `Bo`/`Rs`.

        Raises
        ------
        InputValidationError
            Если давления <= 0 или температура ниже абсолютного нуля.
        """
        validate_positive_pressure(pressure_arr_bar, name='pressure_arr_bar')
        validate_positive_pressure(reservoir_pressure_bar, name='reservoir_pressure_bar')
        validate_temperature_celsius(reservoir_temperature_c, name='reservoir_temperature_c')
        calc = DLE(self._composition_copy(), list(pressure_arr_bar),
                   reservoir_pressure_bar, reservoir_temperature_c)
        calc.calculate()
        df = calc._dle_df

        self._model.result_store_object.add(
            module='Experiments.dle',
            params={
                'pressure_arr_bar': list(pressure_arr_bar),
                'reservoir_pressure_bar': reservoir_pressure_bar,
                'reservoir_temperature_c': reservoir_temperature_c,
            },
            data=df,
        )
        return df

    def cce(self, pressure_arr_bar: list, reservoir_temperature: float) -> pd.DataFrame:
        """
        Дифференциация при постоянном составе (см. `calc_core/Experiments/CCE.py`).

        Parameters
        ----------
        pressure_arr_bar : list
            Ступени давления, бар.
        reservoir_temperature : float
            Температура эксперимента, K (не °C — так требует сам `CCE.__init__`).

        Returns
        -------
        pandas.DataFrame
            Таблица по ступеням давления с колонками `v/v_res`, `v/v_sat`, `Compressibility`.

        Raises
        ------
        InputValidationError
            Если давления <= 0 или температура <= 0 K.
        """
        validate_positive_pressure(pressure_arr_bar, name='pressure_arr_bar')
        validate_temperature_kelvin(reservoir_temperature, name='reservoir_temperature')
        calc = CCE(self._composition_copy(), list(pressure_arr_bar), reservoir_temperature)
        df = calc.calculate()

        self._model.result_store_object.add(
            module='Experiments.cce',
            params={
                'pressure_arr_bar': list(pressure_arr_bar),
                'reservoir_temperature': reservoir_temperature,
            },
            data=df,
        )
        return df

    def separator(self, pressure_arr_bar: list, temperature_arr_c: list,
                  reservoir_pressure_bar: float, reservoir_temperature_c: float) -> pd.DataFrame:
        """
        Многоступенчатая сепарация (см. `calc_core/Experiments/SeparatorTest.py`).

        Parameters
        ----------
        pressure_arr_bar : list
            Давления ступеней, бар.
        temperature_arr_c : list
            Температуры ступеней, °C — той же длины, что `pressure_arr_bar`.
        reservoir_pressure_bar : float
            Пластовое давление, бар.
        reservoir_temperature_c : float
            Пластовая температура, °C.

        Returns
        -------
        pandas.DataFrame
            Таблица по ступеням с колонками `Bo`/`Rs`.

        Raises
        ------
        InputValidationError
            Если давления <= 0 или температуры (пластовая/ступеней) ниже абсолютного нуля.
        LengthMismatchError
            Если `pressure_arr_bar` и `temperature_arr_c` разной длины
            (поднимается изнутри `SeparatorTest.calculate()`).
        """
        validate_positive_pressure(pressure_arr_bar, name='pressure_arr_bar')
        validate_positive_pressure(reservoir_pressure_bar, name='reservoir_pressure_bar')
        validate_temperature_celsius(reservoir_temperature_c, name='reservoir_temperature_c')
        for stage_t in temperature_arr_c:
            validate_temperature_celsius(stage_t, name='temperature_arr_c')
        calc = SeparatorTest(self._composition_copy(), list(pressure_arr_bar), list(temperature_arr_c),
                             reservoir_pressure_bar, reservoir_temperature_c)
        calc.calculate()
        df = calc._dle_df

        self._model.result_store_object.add(
            module='Experiments.separator',
            params={
                'pressure_arr_bar': list(pressure_arr_bar),
                'temperature_arr_c': list(temperature_arr_c),
                'reservoir_pressure_bar': reservoir_pressure_bar,
                'reservoir_temperature_c': reservoir_temperature_c,
            },
            data=df,
        )
        return df
