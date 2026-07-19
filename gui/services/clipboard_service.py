"""Parsing helpers for data exchanged with spreadsheet applications."""

import csv
import io
import math
import re


def _header_key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().casefold())


def _number(value: object) -> float | None:
    text = str(value).strip().replace("\u00a0", "")
    if not text:
        return None
    # Excel locale variants: 1,23; 1.234,56; and 1,234.56.
    if "," in text and "." in text:
        text = (text.replace(".", "") if text.rfind(",") > text.rfind(".")
                else text.replace(",", ""))
    elif "," in text:
        text = text.replace(",", ".")
    if text.endswith("%"):
        text = text[:-1].strip()
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def parse_lab_clipboard(
    text: str,
    columns: list[str],
    start_column: int = 0,
) -> tuple[list[list[float | None]], bool]:
    """Parses copied Excel TSV into lab rows and detects an optional header.

    A header may contain the lab column names in any order. Without a header,
    values are assigned positionally starting at ``start_column``; missing
    cells remain blank and extra columns are ignored.  A one-line comma or
    semicolon separated list is accepted as a convenience for clipboard
    providers that flatten a copied spreadsheet column.
    """
    if not text or not columns:
        return [], False
    start_column = max(0, int(start_column))
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    # A normal Excel copy is TSV, but some clipboard bridges turn a copied
    # single column into ``100, 90, 80``.  Do not split a legitimate decimal
    # comma value such as ``1,23``.
    if "\t" not in normalized and "\n" not in normalized:
        one_line = normalized.strip()
        if _number(one_line) is None:
            separator = ";" if ";" in one_line else ","
            if separator in one_line:
                parsed = list(csv.reader(io.StringIO(one_line), delimiter=separator))
                # Treat a flattened single-column payload as vertical data.
                records = [[value] for value in parsed[0]] if len(parsed) == 1 else parsed
            else:
                records = list(csv.reader(io.StringIO(normalized), delimiter="\t"))
        else:
            records = [[one_line]]
    else:
        records = list(csv.reader(io.StringIO(normalized), delimiter="\t"))
    records = [row for row in records if any(str(cell).strip() for cell in row)]
    if not records:
        return [], False

    known = {_header_key(column): index for index, column in enumerate(columns)}
    pressure_index = known.get("pressure")
    if pressure_index is not None:
        known.update({"pressurebar": pressure_index, "pbar": pressure_index,
                      "p": pressure_index})
    header_map: list[int] | None = None
    first = [_header_key(value) for value in records[0]]
    if sum(key in known for key in first) >= 1:
        header_map = [known.get(key, -1) for key in first]
        records = records[1:]

    rows: list[list[float | None]] = []
    for record in records:
        row: list[float | None] = [None] * len(columns)
        for source_index, value in enumerate(record):
            target_index = (
                header_map[source_index]
                if header_map is not None and source_index < len(header_map)
                else start_column + source_index
            )
            if 0 <= target_index < len(columns):
                row[target_index] = _number(value)
        if any(value is not None for value in row):
            rows.append(row)
    return rows, header_map is not None
