from gui.services.clipboard_service import parse_lab_clipboard
from gui.view.table_clipboard import table_to_tsv


def test_parse_lab_clipboard_reorders_header_and_decimal_comma():
    rows, had_header = parse_lab_clipboard(
        "Bo\tPressure, bar\n1,23\t100\n1,10\t90",
        ["pressure", "Bo"],
    )

    assert had_header is True
    assert rows == [[100.0, 1.23], [90.0, 1.10]]


def test_parse_lab_clipboard_supports_positional_single_column():
    rows, had_header = parse_lab_clipboard("100\n90\n", ["pressure", "Bo"])

    assert had_header is False
    assert rows == [[100.0, None], [90.0, None]]


def test_parse_lab_clipboard_accepts_flattened_column_and_target_column():
    rows, had_header = parse_lab_clipboard(
        "100, 90, 80", ["pressure", "Bo"], start_column=1)

    assert had_header is False
    assert rows == [[None, 100.0], [None, 90.0], [None, 80.0]]


def test_parse_lab_clipboard_accepts_flattened_column_without_spaces():
    rows, had_header = parse_lab_clipboard(
        "300,212,5.5,1.6", ["pressure", "Bo"])

    assert had_header is False
    assert rows == [[300.0, None], [212.0, None], [5.5, None], [1.6, None]]


def test_table_to_tsv_is_excel_compatible_and_sanitizes_cells():
    assert table_to_tsv(["A", "B"], [[1, "two\tlines"], [None, 3]]) == (
        "A\tB\n1\ttwo lines\n\t3"
    )
