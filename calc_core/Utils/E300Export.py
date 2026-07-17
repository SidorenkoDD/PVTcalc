"""
Экспорт состава флюида в формат Eclipse 300 (`.inc`).

Портирован из авторского модуля экспорта (старый пакет `calculations.*`) и
адаптирован под текущую модель данных:
- источник — `Composition` (`.composition` — мольные доли, `.composition_data`
  — свойства; BIP лежит в `composition_data['bip']` как вложенный словарь
  `{comp_i: {comp_j: k_ij}}`, а не отдельной `BIPS`-таблицей);
- единицы METRIC — значения пишутся как есть (Pc в бар, Tc в K, MW в г/моль,
  Vc/ω/shift как в БД), без старых конверсий `MW/1000`/`Pc*10`;
- ключ EOS (`MPR`/`SRK`) задаётся снаружи — модели считаются на EOS
  Брусиловского, у которого нет прямого E300-эквивалента, поэтому уравнение
  для дека выбирает пользователь при экспорте.

Пишет ключевые слова EOS, NCOMPS, CNAMES, MW, TCRIT, PCRIT, VCRIT, ACF,
SSHIFT, ZI, BIC. Не импортирует DearPyGui.
"""

import logging

from calc_core.Composition.Composition import Composition

logger = logging.getLogger(__name__)

_EOS_KEYWORDS = ("MPR", "SRK")


class E300Exporter:
    """Сборка E300-дека (`.inc`) из `Composition`."""

    def __init__(self, composition: Composition, eos_keyword: str = "MPR"):
        """
        Parameters
        ----------
        composition : Composition
            Состав с посчитанными `composition_data`.
        eos_keyword : str
            Уравнение состояния для дека: `"MPR"` (modified Peng-Robinson) или
            `"SRK"`. Брусиловский обобщает PR — по умолчанию `MPR`.
        """
        kw = (eos_keyword or "MPR").upper()
        if kw not in _EOS_KEYWORDS:
            raise ValueError(f"EOS для E300 должен быть одним из {_EOS_KEYWORDS}, а не '{eos_keyword}'")
        self._c = composition
        self._eos = kw

    @property
    def _names(self) -> list:
        return list(self._c.composition.keys())

    def _col(self, key: str) -> dict:
        return self._c.composition_data.get(key, {})

    @staticmethod
    def _block(keyword: str, comment: str, lines: list) -> str:
        """Один блок E300: KEYWORD, шапка-комментарий, строки значений, `/`."""
        out = [keyword, "--", f"--{comment}", "--"]
        out.extend(str(x) for x in lines)
        out.append("/")
        return "\n".join(out) + "\n"

    def _scalar_block(self, keyword: str, comment: str, data_key: str,
                      factor: float = 1.0) -> str:
        """Блок «по одному значению свойства `data_key` на компонент»."""
        col = self._col(data_key)
        vals = [col.get(name, 0.0) * factor for name in self._names]
        return self._block(keyword, comment, vals)

    def _bic_block(self) -> str:
        """
        Нижний треугольник матрицы BIP (Eclipse BIC): для компонента i (со
        второго) — коэффициенты с компонентами 0..i-1, каждая строка отдельно.
        """
        names = self._names
        bip = self._col("bip")
        lines = []
        for i in range(1, len(names)):
            row = [bip.get(names[i], {}).get(names[j], 0.0) for j in range(i)]
            lines.append("  ".join(str(v) for v in row))
        return self._block("BIC", "Binary interaction coefficients", lines)

    def export_string(self) -> str:
        """Весь E300-дек как строка."""
        names = self._names
        zi = self._c.composition
        parts = [
            "--E300.inc\n",
            "-- Auto generated file with PVTcalc\n",
            self._block("EOS", "Equation of state", [self._eos]),
            self._block("NCOMPS", "Number of components", [len(names)]),
            self._block("CNAMES", "Component names", names),
            self._scalar_block("MW", "Molecular weight, g/mol", "molar_mass"),
            self._scalar_block("TCRIT", "Critical temperature, K", "critical_temperature"),
            self._scalar_block("PCRIT", "Critical pressure, bar", "critical_pressure"),
            self._scalar_block("VCRIT", "Critical volume", "critical_volume"),
            self._scalar_block("ACF", "Acentric factors", "acentric_factor"),
            self._scalar_block("SSHIFT", "EoS volume shift", "shift_parameter"),
            self._block("ZI", "Overall composition", [zi[name] for name in names]),
            self._bic_block(),
        ]
        return "\n".join(parts)

    def export(self, path: str) -> str:
        """
        Пишет дек в файл. Расширение `.inc` добавляется, если его нет.
        Возвращает фактический путь записанного файла.
        """
        if not path.lower().endswith(".inc"):
            path = path + ".inc"
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.export_string())
        logger.info("E300-дек записан: %s (%d компонентов, EOS=%s)",
                    path, len(self._names), self._eos)
        return path
