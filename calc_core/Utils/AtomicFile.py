"""Атомарная запись файлов, используемых движком и GUI.

Данные сначала пишутся во временный файл в том же каталоге, принудительно
сбрасываются на диск, после чего ``os.replace`` атомарно подменяет целевой
файл. Поэтому авария посередине сериализации не оставляет усечённый JSON.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any


def atomic_write_text(path: str | Path, text: str, *, encoding: str = "utf-8") -> None:
    """Атомарно записывает ``text`` в ``path``."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
                mode="w", encoding=encoding, newline="", delete=False,
                dir=target.parent, prefix=f".{target.name}.", suffix=".tmp") as tmp:
            tmp.write(text)
            tmp.flush()
            os.fsync(tmp.fileno())
            tmp_path = Path(tmp.name)
        os.replace(tmp_path, target)
        tmp_path = None
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)


def atomic_write_json(path: str | Path, data: Any, *, default=None) -> None:
    """Сериализует ``data`` и атомарно записывает JSON в ``path``."""
    text = json.dumps(data, indent=2, default=default, ensure_ascii=False)
    atomic_write_text(path, text + "\n")
