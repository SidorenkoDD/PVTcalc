"""Golden-master на PhaseEnvelopeSSM — каноническую реализацию фазовой огибающей.

Зачем этот файл появился (2026-07-19): `calc_core/PhaseEnvelope/` — область с
наибольшим числом конкурирующих реализаций и при этом 10 из 11 модулей не имели
ни одного прямого теста. SSM при этом основная реализация и её вот-вот будут
править (протаскивание отмены/прогресса, перевод параллельного прохода с
`joblib.Parallel` на явный futures-цикл). Править непокрытый численный модуль —
верный способ молча сломать физику.

Как и `test_flash.py`, это НЕ проверка физической правильности (независимого
эталона нет), а фиксация текущего вывода движка.

Что здесь важно понимать про воспроизводимость — замерено, а не предположено:

* `calculate()` полностью детерминирован при фиксированной сетке T. Но значения
  зависят от САМОЙ сетки: continuation передаёт давление с предыдущего шага как
  затравку, поэтому при T=100 °C шаг 10 и шаг 20 дают чуть разные числа
  (178.8374 против 178.8552). Эталон закреплён под конкретную сетку.
* `calculate_parallel(n_jobs=1)` совпадает с `calculate()` до машинной точности
  (~4e-16): один воркер = один кусок = тот же continuation. Это самый сильный
  инвариант в файле и, в отличие от многоворкерного, он не зависит от числа
  ядер машины.
* `calculate_parallel(n_jobs>1)` расходится с последовательным на ~1e-4
  относительных. Это не ошибка: куски разрывают continuation, а солвер и так
  бисектится с `bracket_tol_rel=1e-3`. Поэтому абсолютные значения для
  многоворкерного прохода не закрепляются — на CI (2-4 ядра) разбиение другое,
  чем на машине автора (20).

Входной состав берётся из версионированного `tests/fixtures/models.json` и не
меняется при работе пользователя с корневой базой. Это regression baseline
текущего движка, а не независимый физический эталон.
"""

import numpy as np
import pytest

from calc_core.PhaseEnvelope.PhaseEnvelopeSuccessiveSubstitution import PhaseEnvelopeSSM

# Сетка подобрана под CI: ~0.6 c последовательно, ~0.7 c параллельно.
T_MIN_C, T_MAX_C, T_STEP_C, P_MAX_BAR = 60.0, 140.0, 20.0, 400.0

# calculate() на KRSNL_PVTSIM, снято 2026-07-19.
EXPECTED_BUBBLE_BAR = {
    60.0: 149.04357132419398,
    80.0: 164.89400581365170,
    100.0: 178.85524634494425,
    120.0: 190.74981888018905,
    140.0: 200.64124015219886,
}

# Допуск сравнения последовательного с многоворкерным: порядка собственного
# bracket_tol_rel солвера (1e-3), а не машинной точности.
SOLVER_REL_TOL = 1e-3


def _ssm(composition):
    return PhaseEnvelopeSSM(composition, T_MIN_C, T_MAX_C, T_STEP_C, P_MAX_BAR)


def test_ssm_sequential_golden_master(krsnl_composition):
    """Точные значения последовательного марша по фиксированной сетке T."""
    df = _ssm(krsnl_composition).calculate()

    assert list(df.columns) == ["Temp_C", "Bubble_bar", "Dew_upper_bar", "Dew_lower_bar"]
    assert df["Temp_C"].tolist() == list(EXPECTED_BUBBLE_BAR)

    for temp, expected in EXPECTED_BUBBLE_BAR.items():
        got = float(df.loc[df["Temp_C"] == temp, "Bubble_bar"].iloc[0])
        assert got == pytest.approx(expected, rel=1e-9), f"T={temp} °C"

    # У этого состава в этом диапазоне SSM находит только ветку bubble —
    # ретроградная область не захвачена (замечание автора, docs/GUI.md СТАТУС-4).
    assert df["Dew_upper_bar"].isna().all()
    assert df["Dew_lower_bar"].isna().all()


def test_ssm_sequential_is_deterministic(krsnl_composition):
    """Повторный прогон на том же составе даёт идентичный результат."""
    first = _ssm(krsnl_composition).calculate()
    second = _ssm(krsnl_composition).calculate()

    assert first.equals(second)


def test_ssm_parallel_single_worker_matches_sequential(krsnl_composition):
    """n_jobs=1 — один кусок, тот же continuation, совпадение до машинной точности.

    Самый сильный инвариант файла: не зависит от числа ядер, поэтому одинаково
    строг на машине автора и на CI. Именно он ловит подмену логики сборки
    результатов при переходе на futures-цикл.
    """
    sequential = _ssm(krsnl_composition).calculate()
    parallel = _ssm(krsnl_composition).calculate_parallel(n_jobs=1)

    assert parallel["Temp_C"].tolist() == sequential["Temp_C"].tolist()
    np.testing.assert_allclose(
        parallel["Bubble_bar"].to_numpy(dtype=float),
        sequential["Bubble_bar"].to_numpy(dtype=float),
        rtol=1e-12,
    )


@pytest.mark.parametrize("n_jobs", [2, 3, -1])
def test_ssm_parallel_preserves_temperature_order(krsnl_composition, n_jobs):
    """Температуры идут строго по возрастанию и совпадают с заданной сеткой.

    Ловушка, ради которой тест написан: `joblib.Parallel` отдаёт результаты в
    порядке подачи задач, а `concurrent.futures.as_completed` — в порядке
    завершения. При переводе параллельного прохода на futures куски, собранные
    без явной сортировки, дадут перемешанную огибающую — численно правдоподобную
    и потому незаметную глазом на графике.
    """
    df = _ssm(krsnl_composition).calculate_parallel(n_jobs=n_jobs)

    temps = df["Temp_C"].to_numpy(dtype=float)
    assert temps.tolist() == list(EXPECTED_BUBBLE_BAR)
    assert np.all(np.diff(temps) > 0), f"температуры не по возрастанию: {temps.tolist()}"


def test_ssm_parallel_matches_sequential_within_solver_tolerance(krsnl_composition):
    """Многоворкерный проход сходится с последовательным в пределах допуска солвера.

    Разрыв continuation на границах кусков сдвигает результат бисекции примерно
    на 1e-4 относительных — меньше собственного `bracket_tol_rel=1e-3` солвера.
    Проверяем именно порядок величины расхождения: рост на порядки означал бы,
    что разбиение на куски стало ломать сходимость.
    """
    sequential = _ssm(krsnl_composition).calculate()
    parallel = _ssm(krsnl_composition).calculate_parallel(n_jobs=-1)

    np.testing.assert_allclose(
        parallel["Bubble_bar"].to_numpy(dtype=float),
        sequential["Bubble_bar"].to_numpy(dtype=float),
        rtol=SOLVER_REL_TOL,
    )
