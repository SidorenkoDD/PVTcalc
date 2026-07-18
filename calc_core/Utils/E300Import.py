"""
Импорт состава флюида из формата Eclipse 300 (`.inc`).

Обратная операция к `E300Export`: разбирает дек по ключевым словам
(CNAMES/ZI/MW/TCRIT/PCRIT/VCRIT/ACF/SSHIFT/BIC/EOS) в набор компонентов с их
свойствами. Единицы — METRIC (как пишет `E300Exporter`): Pc в бар, Tc в K,
MW в г/моль; конверсий нет.

Парсер устойчив к комментариям (`--…`) и лишним пробелам: блок начинается
строкой-ключевым словом и заканчивается строкой `/`; строки-комментарии
внутри блока игнорируются, остальные — значения (в BIC на строке несколько
значений — нижний треугольник матрицы BIP).

Не строит `Composition` сам (это делает сервисный слой GUI) — только парсит.
"""

import logging
import math

logger = logging.getLogger(__name__)

# Ключевое слово E300 -> имя свойства в composition_data.
_PROP_KEYWORDS = {
    "MW": "molar_mass",
    "TCRIT": "critical_temperature",
    "PCRIT": "critical_pressure",
    "VCRIT": "critical_volume",
    "ACF": "acentric_factor",
    "SSHIFT": "shift_parameter",
}
_ALL_KEYWORDS = set(_PROP_KEYWORDS) | {"EOS", "NCOMPS", "CNAMES", "ZI", "BIC"}


def _read_blocks(text: str) -> dict:
    """
    Разбирает текст дека в `{keyword: [строки-значений]}`.

    Строки-комментарии (`--…`) и пустые пропускаются; блок идёт от строки с
    ключевым словом до строки `/`.
    """
    blocks: dict[str, list] = {}
    current = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("--"):
            continue
        if line == "/":
            current = None
            continue
        if current is None:
            # ожидаем ключевое слово (первый токен)
            kw = line.split()[0].upper()
            if kw in _ALL_KEYWORDS:
                current = kw
                blocks[current] = []
            # незнакомые ключевые слова просто пропускаем до их «/»
            else:
                current = "_SKIP_"
                blocks.setdefault("_SKIP_", [])
            continue
        blocks[current].append(line)
    blocks.pop("_SKIP_", None)
    return blocks


def _floats(lines: list) -> list:
    """Плоский список float из строк блока (в строке может быть несколько)."""
    out = []
    for ln in lines:
        for tok in ln.replace(",", " ").split():
            out.append(float(tok))
    return out


def parse_e300(path: str) -> dict:
    """
    Разбирает E300-дек в структуру для восстановления состава.

    Returns
    -------
    dict
        ``{"eos": str|None, "names": [str], "zi": {name: float},
        "props": {data_key: {name: float}}, "bip": {ci: {cj: float}}}``.
        `props` содержит только присутствующие в деке свойства.

    Raises
    ------
    ValueError
        Если нет обязательных блоков CNAMES/ZI или длины не согласованы.
    """
    with open(path, "r", encoding="utf-8") as f:
        blocks = _read_blocks(f.read())

    if "CNAMES" not in blocks or "ZI" not in blocks:
        raise ValueError("E300 deck must contain CNAMES and ZI blocks")

    warnings: list[str] = []
    names = [ln.split()[0] for ln in blocks["CNAMES"]]
    if not names or len(set(names)) != len(names):
        raise ValueError("CNAMES must contain unique component names")
    if blocks.get("NCOMPS"):
        ncomps = int(_floats(blocks["NCOMPS"])[0])
        if ncomps != len(names):
            raise ValueError(f"NCOMPS ({ncomps}) does not match CNAMES ({len(names)})")
    zi_vals = _floats(blocks["ZI"])
    if len(zi_vals) != len(names):
        raise ValueError(
            f"ZI count ({len(zi_vals)}) does not match CNAMES ({len(names)})")
    if any(not math.isfinite(v) or v < 0.0 for v in zi_vals) or sum(zi_vals) <= 0.0:
        raise ValueError("ZI values must be finite, non-negative and have a positive sum")
    zi = {n: v for n, v in zip(names, zi_vals)}

    props: dict[str, dict] = {}
    for kw, data_key in _PROP_KEYWORDS.items():
        if kw not in blocks:
            continue
        vals = _floats(blocks[kw])
        if len(vals) != len(names):
            logger.warning("E300 %s: %d значений при %d компонентах — пропущено",
                           kw, len(vals), len(names))
            warnings.append(f"{kw}: {len(vals)} values for {len(names)} components; block skipped")
            continue
        if any(not math.isfinite(v) for v in vals):
            warnings.append(f"{kw}: non-finite values; block skipped")
            continue
        props[data_key] = {n: v for n, v in zip(names, vals)}

    # BIC -> симметричная матрица BIP (нижний треугольник: строка i -> j<i)
    bip = {ci: {cj: 0.0 for cj in names} for ci in names}
    if "BIC" in blocks:
        if len(blocks["BIC"]) < max(0, len(names) - 1):
            warnings.append("BIC: incomplete lower-triangle matrix; missing values set to 0")
        for i, ln in enumerate(blocks["BIC"], start=1):
            if i >= len(names):
                break
            row = [float(t) for t in ln.replace(",", " ").split()]
            if len(row) != i:
                warnings.append(f"BIC row {i + 1}: expected {i} values, got {len(row)}")
            for j, val in enumerate(row):
                if j < i and math.isfinite(val):
                    bip[names[i]][names[j]] = val
                    bip[names[j]][names[i]] = val

    eos = blocks["EOS"][0].split()[0].upper() if blocks.get("EOS") else None
    logger.info("E300-дек разобран: %d компонентов, EOS=%s, свойств=%d",
                len(names), eos, len(props))
    return {"eos": eos, "names": names, "zi": zi, "props": props,
            "bip": bip, "warnings": warnings}
