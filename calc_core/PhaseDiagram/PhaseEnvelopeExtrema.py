"""
Криконденбара и крикондентерма — экстремумы фазовой огибающей.

Обе величины ищутся в два шага: (1) грубое приближение — argmax по уже
построенной таблице `PhaseEnvelopeSSM` (`Temp_C`/`Bubble_bar`/`Dew_upper_bar`),
(2) уточнение — 1D-максимизация вдоль ВЕРХНЕЙ (внешней) ветки огибающей в
узкой скобке вокруг грубого приближения (см. `_upper_branch`). Обе величины
лежат именно на этой ветке, не на `Dew_lower_bar` (внутренняя ретроградная
петля — её T/P-диапазон является подмножеством внешней ветки, см. докстринг
`PhaseEnvelopeSuccessiveSubstitution.py`).

Криконденбара (максимум P) — касательная к огибающей там горизонтальна
(dP/dT=0), поэтому это прямая максимизация P(T); существующий T->P решатель
(`SaturationPointSSM.find_bubble_point`) подходит без изменений — тот же
решатель, что уже используется в `PhaseEnvelopeSSM`.

Крикондентерма (максимум T) — касательная вертикальна (dT/dP=0). Максимизация
P(T) вблизи такой точки численно плохо обусловлена (маленькое изменение T ->
большое изменение P), поэтому уточнение использует T->P-обратный решатель
(`SaturationTemperatureSSM.find_dew_temperature`, P задано, T ищется) — см. её
докстринг за подробным обоснованием.
"""

import logging

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from calc_core.Composition.Composition import Composition
from calc_core.PhaseDiagram.PhaseEnvelopeSuccessiveSubstitution import (
    PhaseEnvelopeSSM,
    SaturationPointSSM,
    SaturationTemperatureSSM,
)

logger = logging.getLogger(__name__)


def _upper_branch(envelope_df: pd.DataFrame) -> pd.DataFrame:
    """
    Объединённая "верхняя" (внешняя) ветка огибающей — `Bubble_bar` и
    `Dew_upper_bar` склеены через `combine_first` (тот же приём, что уже
    используется в `PhaseEnvelopeSSM.plot()`): это одна физическая кривая,
    разбитая на два столбца только классификацией Михельсена по разные
    стороны критической точки состава.

    Returns
    -------
    pd.DataFrame
        Столбцы `Temp_C`/`P_bar`, отсортировано по `Temp_C`, без NaN.
    """
    df = envelope_df.sort_values('Temp_C').reset_index(drop=True)
    p = df['Bubble_bar'].combine_first(df['Dew_upper_bar'])
    mask = p.notna()
    return pd.DataFrame({
        'Temp_C': df.loc[mask, 'Temp_C'].to_numpy(),
        'P_bar': p[mask].to_numpy(),
    })


def _build_coarse_envelope(
    composition: Composition, p_max_bar: float, p_min_bar: float,
    t_min_c: float, t_max_c: float | None, t_step_c: float,
) -> pd.DataFrame:
    if t_max_c is None:
        raise ValueError(
            "Нужно передать либо envelope_df (уже построенную огибающую), "
            "либо t_max_c для построения грубой огибающей внутри."
        )
    logger.info("envelope_df не передан — строим грубую огибающую (SSM) внутри")
    return PhaseEnvelopeSSM(
        composition, t_min_c, t_max_c, t_step_c, p_max_bar, p_min_bar,
    ).calculate()


class CricondenbarCalculator:
    """
    Криконденбара — точка максимального давления на фазовой огибающей.
    См. докстринг модуля за обоснованием алгоритма.
    """

    def __init__(
        self,
        composition: Composition,
        p_max_bar: float,
        envelope_df: pd.DataFrame | None = None,
        p_min_bar: float = 0.1,
        t_min_c: float = 0.0,
        t_max_c: float | None = None,
        t_step_c: float = 10.0,
    ):
        self.composition = composition
        self.p_max_bar = p_max_bar
        self.p_min_bar = p_min_bar

        if envelope_df is None:
            envelope_df = _build_coarse_envelope(
                composition, p_max_bar, p_min_bar, t_min_c, t_max_c, t_step_c,
            )

        self.envelope_df = envelope_df
        self.result = None

    def _p_at_t(self, t_c: float, p_guess: float) -> float | None:
        t_k = t_c + 273.15
        point = SaturationPointSSM(
            self.composition, t_k, self.p_max_bar, self.p_min_bar,
            p_guess_bubble=p_guess,
        ).find_bubble_point()
        return point.pressure if point is not None else None

    def calculate(self) -> dict:
        """
        Returns
        -------
        dict
            `{'T': K, 'T_C': °C, 'P': бар}`.
        """
        upper = _upper_branch(self.envelope_df)
        if upper.empty:
            raise RuntimeError("Огибающая пуста — не на чем искать криконденбару.")

        p_arr = upper['P_bar'].to_numpy()
        t_arr = upper['Temp_C'].to_numpy()
        idx_coarse = int(np.argmax(p_arr))
        t_coarse = float(t_arr[idx_coarse])
        p_coarse = float(p_arr[idx_coarse])

        unique_t = np.sort(np.unique(t_arr))
        t_step = float(np.median(np.diff(unique_t))) if unique_t.size > 1 else 5.0

        t_lo = max(float(t_arr.min()), t_coarse - t_step)
        t_hi = min(float(t_arr.max()), t_coarse + t_step)

        def neg_p(t_c: float) -> float:
            p = self._p_at_t(t_c, p_guess=p_coarse)
            return -p if p is not None else 1e6

        t_refined, p_refined = t_coarse, p_coarse
        if t_hi > t_lo:
            res = minimize_scalar(neg_p, bounds=(t_lo, t_hi), method='bounded', options={'xatol': 1e-3})
            if res.success:
                p_at_res = self._p_at_t(float(res.x), p_guess=p_coarse)
                if p_at_res is not None:
                    t_refined, p_refined = float(res.x), p_at_res

        self.result = {'T': t_refined + 273.15, 'T_C': t_refined, 'P': p_refined}
        logger.info("Криконденбара: T=%.2f °C, P=%.4f бар", t_refined, p_refined)
        return self.result


class CricondenthermCalculator:
    """
    Крикондентерма — точка максимальной температуры на фазовой огибающей.
    См. докстринг модуля за обоснованием алгоритма.
    """

    def __init__(
        self,
        composition: Composition,
        p_max_bar: float,
        envelope_df: pd.DataFrame | None = None,
        p_min_bar: float = 0.1,
        t_min_c: float = 0.0,
        t_max_c: float | None = None,
        t_step_c: float = 10.0,
    ):
        self.composition = composition
        self.p_max_bar = p_max_bar
        self.p_min_bar = p_min_bar

        if envelope_df is None:
            envelope_df = _build_coarse_envelope(
                composition, p_max_bar, p_min_bar, t_min_c, t_max_c, t_step_c,
            )

        self.envelope_df = envelope_df
        self.result = None

    def _t_at_p(self, p_bar: float, t_guess_k: float, t_min_k: float, t_max_k: float) -> float | None:
        point = SaturationTemperatureSSM(
            self.composition, p_bar, t_max_k, t_min_k,
            t_guess_dew=t_guess_k,
        ).find_dew_temperature(t_lower_bound=t_min_k)
        return point.temperature if point is not None else None

    def calculate(self) -> dict:
        """
        Returns
        -------
        dict
            `{'T': K, 'T_C': °C, 'P': бар}`.
        """
        upper = _upper_branch(self.envelope_df)
        if upper.empty:
            raise RuntimeError("Огибающая пуста — не на чем искать крикондентерму.")

        p_arr = upper['P_bar'].to_numpy()
        t_arr = upper['Temp_C'].to_numpy()
        idx_coarse = int(np.argmax(t_arr))
        t_coarse = float(t_arr[idx_coarse])
        p_coarse = float(p_arr[idx_coarse])

        unique_p = np.sort(np.unique(p_arr))
        p_step = float(np.median(np.diff(unique_p))) if unique_p.size > 1 else max(1.0, 0.05 * p_coarse)

        p_lo = max(self.p_min_bar, p_coarse - p_step)
        p_hi = min(self.p_max_bar, p_coarse + p_step)

        # Безопасный диапазон T для решателя (см. SaturationTemperatureSSM) —
        # шире, чем узкая скобка guess'а (30%, см. find_dew_temperature), но
        # сам по себе не участвует в scipy-скобке (та — по P).
        t_guess_k = t_coarse + 273.15
        t_min_k = max(1.0, t_guess_k * 0.5)
        t_max_k = t_guess_k * 1.5

        def neg_t(p_bar: float) -> float:
            t = self._t_at_p(p_bar, t_guess_k=t_guess_k, t_min_k=t_min_k, t_max_k=t_max_k)
            return -t if t is not None else 1e6

        t_refined_k, p_refined = t_guess_k, p_coarse
        if p_hi > p_lo:
            res = minimize_scalar(neg_t, bounds=(p_lo, p_hi), method='bounded', options={'xatol': 1e-3})
            if res.success:
                t_at_res = self._t_at_p(float(res.x), t_guess_k=t_guess_k, t_min_k=t_min_k, t_max_k=t_max_k)
                if t_at_res is not None:
                    t_refined_k, p_refined = t_at_res, float(res.x)

        t_refined_c = t_refined_k - 273.15
        self.result = {'T': t_refined_k, 'T_C': t_refined_c, 'P': p_refined}
        logger.info("Крикондентерма: T=%.2f °C, P=%.4f бар", t_refined_c, p_refined)
        return self.result
