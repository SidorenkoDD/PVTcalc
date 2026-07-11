"""
Загрузчик статической БД компонентов (`Utils/DB.json`).

Устойчив к тому, откуда запущен код (ноутбук, pytest, установленный пакет) —
ищет файл по нескольким кандидатным путям (см. `_find_db_file`), а не по
жёстко заданному относительному пути. Используется в `Composition.py` как
модульная переменная `db = JsonDBReader().load_database('DB.json')`,
загружается один раз при импорте.
"""

import json
import os
import sys
from pathlib import Path
from calc_core.Utils.BaseClasses import Reader


class JsonDBReader(Reader):
    """Читает JSON-файл базы данных компонентов, разрешая путь к нему автоматически."""

    def _find_db_file(self, db_filename):
        """
        Поиск DB.json с учетом различных сценариев запуска.

        Перебирает кандидатные пути в порядке: рядом с этим файлом
        (`calc_core/Utils/`) → на уровень выше → текущая рабочая директория
        → рядом с исполняемым скриптом (`sys.argv[0]`) → переменная
        окружения `DB_PATH` → до 4 уровней вверх по дереву директорий.
        Возвращает первый существующий путь.

        Именно первый кандидат ("рядом с этим файлом") и делает загрузку
        рабочей после `pip install`/`pip install -e .` — `DB.json` бандлится
        в пакет рядом с исходниками (`package-data` в `pyproject.toml`), и
        именно там он и находится, независимо от cwd.

        Parameters
        ----------
        db_filename : str
            Имя файла (обычно `'DB.json'`).

        Returns
        -------
        Path

        Raises
        ------
        FileNotFoundError
            Если файл не найден ни по одному кандидатному пути.
        """
        current_file_path = Path(__file__).absolute()
        
        possible_paths = [
            # Рядом с текущим файлом
            current_file_path.parent / db_filename,
            # В корне проекта (предполагаем, что текущий файл в подпапке)
            current_file_path.parent.parent / db_filename,
            # В рабочей директории
            Path.cwd() / db_filename,
            # Рядом с исполняемым скриптом
            Path(sys.argv[0]).parent / db_filename,
            # Из переменной окружения
            Path(os.getenv('DB_PATH')) if os.getenv('DB_PATH') else None,
        ]
        
        # Добавляем поиск вверх по иерархии директорий
        for level in range(1, 5):  # Проверяем до 4 уровней вверх
            parent_dir = current_file_path.parents[level]
            possible_paths.append(parent_dir / db_filename)
        
        # Убираем None и дубликаты
        possible_paths = [path for path in possible_paths if path is not None]
        possible_paths = list(dict.fromkeys(possible_paths))  # Remove duplicates
        
        for path in possible_paths:
            if path.exists():
                #print(f"Found DB.json at: {path}")
                return path
        
        
        raise FileNotFoundError(f"{db_filename} not found in any expected location")
    
    def load_database(self, db_filename: str = 'DB.json'):
        """
        Находит и загружает JSON-файл БД.

        Parameters
        ----------
        db_filename : str, optional
            Имя файла. По умолчанию `'DB.json'` (сам по себе неиспользуемый
            дефолт — живой код всегда вызывает с явным `'DB.json'`).

        Returns
        -------
        dict
            Распарсенный JSON (структура — см. `Utils/DB.json`:
            `available_components`, `sequence_number`, `molar_mass`, `Tb`,
            `critical_pressure`/`critical_temperature`/`acentric_factor`/
            `critical_volume`, `bip_pr`/`bip_srk`, `carbon_flag`,
            `c5_plus_flag`/`c7_plus_flag`, `Kw` и т.д.).
        """
        db_path = self._find_db_file(db_filename)
        with open(db_path, 'r', encoding='utf-8') as f:
            return json.load(f)
