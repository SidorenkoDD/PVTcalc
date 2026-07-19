"""Small helpers for copying read-only GUI tables as Excel-compatible TSV."""

from collections.abc import Iterable


def _cell_text(value: object) -> str:
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def table_to_tsv(columns: Iterable[object], rows: Iterable[Iterable[object]]) -> str:
    """Serializes a rectangular table to tab-separated text for Excel/Sheets."""
    lines = ["\t".join(_cell_text(value) for value in columns)]
    lines.extend("\t".join(_cell_text(value) for value in row) for row in rows)
    return "\n".join(lines)
