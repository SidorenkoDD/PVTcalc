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


def parse_row_selection(value: str, total_rows: int) -> list[int]:
    """Parses a one-based row list such as ``1, 3-5`` in table order."""
    if total_rows < 1:
        return []
    text = value.strip()
    if not text:
        return list(range(total_rows))
    selected: set[int] = set()
    for chunk in text.split(","):
        part = chunk.strip()
        if not part:
            continue
        left, separator, right = part.partition("-")
        try:
            first = int(left.strip())
            last = int(right.strip()) if separator else first
        except ValueError as exc:
            raise ValueError("Use row numbers such as 1, 3-5") from exc
        if not separator:
            last = first
        if first < 1 or last < 1 or first > last or last > total_rows:
            raise ValueError(f"Rows must be between 1 and {total_rows}")
        selected.update(range(first - 1, last))
    if not selected:
        raise ValueError("Select at least one row")
    return sorted(selected)


def selected_table(columns: Iterable[object], rows: Iterable[Iterable[object]],
                   row_indices: Iterable[int], column_indices: Iterable[int],
                   ) -> tuple[list[object], list[list[object]]]:
    """Returns a rectangular subset while preserving the original table order."""
    labels = list(columns)
    values = [list(row) for row in rows]
    selected_columns = [index for index in column_indices
                        if 0 <= index < len(labels)]
    if not selected_columns:
        raise ValueError("Select at least one column")
    selected_rows = [index for index in row_indices if 0 <= index < len(values)]
    return (
        [labels[index] for index in selected_columns],
        [[values[row][column] if column < len(values[row]) else ""
          for column in selected_columns] for row in selected_rows],
    )
