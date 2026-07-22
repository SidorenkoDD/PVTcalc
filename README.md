# PVTcalc

PVTcalc — Python-движок компонентного PVT-моделирования и desktop-интерфейс
для работы с моделями нефтяных, газовых и газоконденсатных флюидов.

Проект поддерживает полный рабочий путь: состав флюида → свойства компонентов
и C7+ → EOS Брусиловского → тест стабильности → одно-/двухфазный Flash →
CCE/DLE/Separator → фазовая огибающая P–T. Интерфейс построен на DearPyGui и
работает непосредственно с Python-ядром без отдельного сервера.

> Проект находится в активной разработке. Flash, эксперименты и GUI уже
> используются, но критическая точка, единый Psat-контракт и долговременное
> хранение результатов ещё дорабатываются. Ограничения перечислены ниже.

## Что уже работает

Расчётное ядро:

- состав и свойства компонентов, BIP и корреляции C7+;
- кубическое EOS Брусиловского;
- тест двухфазной стабильности и идентификация однофазного состояния;
- двухфазный Flash с Rachford–Rice и Ньютоном по фугитивностям;
- CCE, DLE и многоступенчатая сепарация;
- фазовая огибающая методами SSM и Grid;
- импорт и экспорт Eclipse 300;
- проверяемое и атомарно обновляемое хранилище моделей.

GUI:

- проекты с несколькими независимыми моделями-вариантами;
- создание состава вручную, импорт из Excel и E300;
- редактирование состава, свойств компонентов, корреляций, EOS и BIP;
- история Flash, CCE, DLE, Separator и Phase envelope;
- фоновые расчёты, статусы `EMPTY/FRESH/STALE/RUNNING` и ленивая инвалидация;
- таблицы, графики, лабораторные точки и копирование данных в Excel;
- undo/redo, сохранение модели и восстановление workspace-сессии.

## Требования

- Python **3.11 или новее**;
- Windows — основная среда desktop GUI;
- Linux поддерживается тестами, GUI-тесты в CI запускаются через Xvfb.

Основные зависимости устанавливаются из `pyproject.toml`. Имя дистрибутива —
`pvtcalc`, импортируемые пакеты — `calc_core` и `gui`.

## Быстрый старт GUI

Из корня репозитория в PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[gui]"
python -m gui
```

По умолчанию приложение читает `models.json` и сохраняет `gui_session.json` в
текущем рабочем каталоге. Другую базу моделей можно передать явно:

```powershell
python -m gui --db "D:\PVT\models.json" --log-level INFO
```

Доступные параметры:

```powershell
python -m gui --help
```

### Первый рабочий сценарий

1. На экране **Projects** создайте модель вручную либо импортируйте Excel/E300.
2. Откройте проект двойным кликом или клавишей Enter.
3. Проверьте **Composition**: мольные доли, свойства, корреляции и BIP.
4. Создайте Flash или эксперимент в левом дереве и задайте условия.
5. После изменения состава узлы с прежними результатами станут `STALE` — их
   нужно пересчитать осознанно.
6. Сохраните модель через **File → Save model** или `Ctrl+S`.

## Установка для разработки

```powershell
python -m pip install -e ".[dev]"
```

Опционально для ноутбуков:

```powershell
python -m pip install -e ".[notebook]"
```

Editable-установка нужна, чтобы `calc_core` и `gui` импортировались независимо
от текущей директории, а изменения исходников подхватывались без переустановки.

## Использование расчётного ядра

Для одно- и двухфазного Flash используйте единый фасад
`CompositionalModel.flash()`. Давление передаётся в барах, температура — в °C:

```python
from calc_core.Composition.Composition import Composition
from calc_core.CompositionalModel.CompositionalModel import CompositionalModel

database = Composition.from_db("models.json")
print(database.list_models())

saved = getattr(database, database.list_models()[0])
composition = saved.new_composition(saved.composition, deep_copy=True)
model = CompositionalModel(composition)

result = model.flash(p_bar=100.0, t_celsius=110.0)
print(result.is_two_phase, result.phase_type)
print(result.vapor.mole_fraction, result.liquid.mole_fraction)
```

Фасад сам выполняет глубокую копию переданного состава: низкоуровневый Flash
меняет температуру и пересчитывает температурно-зависимые параметры рабочей
копии, но исходная модель остаётся неизменной.

## Данные и безопасность записи

`models.json` — общая база проектов и моделей. Живые пути чтения/записи проходят
через `calc_core/Utils/ModelStore.py`, который обеспечивает:

- валидацию схемы и миграцию старых записей в памяти;
- межпроцессный lock и optimistic concurrency;
- атомарную замену файла;
- три ротационные резервные копии: `.bak`, `.bak.1`, `.bak.2`;
- отказ от автоматической перезаписи повреждённой базы.

`gui_session.json` хранит только интерфейс: последний выбор, размер окна и
однократно ещё не мигрированные legacy-workspace. Состав, параметры запусков,
лабораторные точки, результаты и их provenance сохраняются вместе с моделью в
`models.json` после **Save model** или выбора **Save all** при переходе в
Projects / закрытии через **Project → Exit application...**. Таблица Projects
показывает только сохранённую картину и время последнего сохранения.

Тесты не читают пользовательский `models.json`: их входные модели закреплены в
`tests/fixtures/models.json`, поэтому работа в GUI не меняет regression baseline.

## Архитектура

```text
gui/view/                    DearPyGui, экраны и отрисовка
        ↓ commands/events
gui/app_state.py             проекты, модели, вычислительный граф, undo/redo
gui/calculation_coordinator  жизненный цикл фоновой задачи
gui/services/                валидация и адаптеры GUI ↔ calc_core
        ↓
calc_core/CompositionalModel фасады экспериментов и фазовой огибающей
calc_core/VLE                Flash и фазовое равновесие
calc_core/PhaseStability     стабильность и тип одной фазы
calc_core/EOS                EOS Брусиловского
calc_core/Composition        состав, свойства и BIP
```

Интерактивная подробная карта находится в
[`docs/diagrams/module-interactions.html`](docs/diagrams/module-interactions.html).

## Проверки

Полный тестовый набор:

```powershell
python -m pytest -p no:cacheprovider -q
```

Статический анализ:

```powershell
python -m ruff check gui tests
python -m mypy gui
```

CI выполняет pytest на Python 3.11, 3.12 и 3.13, а также отдельно запускает
Ruff и mypy. Численные regression-тесты фиксируют текущий результат движка;
они не считаются внешней верификацией. Официальная сверка с PVTSim для состава
KRSNLN будет добавлена после получения чистых исходных данных.

## Текущие ограничения

- критическая точка недостоверна для многокомпонентных составов и поэтому не
  показывается в GUI;
- GUI-огибающая использует SSM/Grid, тогда как CCE/DLE/Separator пока получают
  Psat через переходный `new_methodv2.SaturationPressure`;
- CCE заметно тяжелее и менее устойчива, чем DLE/Separator;
- Cancel кооперативно останавливает длинные циклы в безопасных точках; старые
  атомарные участки могут завершить текущий внутренний шаг;
- константы и допуски в Settings доступны только для чтения: `EngineConfig`
  уже используется основным Flash-путём, но редактируемый профиль модели ещё
  не подключён;
- CVD, лампинг/делампинг, EOS-regression и mixing-rule заморожены до отдельного
  решения автора.

## Документация и план

- [`docs/ROADMAP.md`](docs/ROADMAP.md) — актуальный порядок будущей работы,
  масштабы задач, зависимости и критерии готовности;
- [`docs/BACKLOG.md`](docs/BACKLOG.md) — история задач и выполненных решений;
- [`docs/NUMERICAL_VERIFICATION.md`](docs/NUMERICAL_VERIFICATION.md) —
  канонические методы фазовой огибающей, допуски, статусы NaN и fallback;
- [`docs/PLUS_COMPONENTS.md`](docs/PLUS_COMPONENTS.md) — каталог корреляций
  C7+, входные/выходные единицы, исполняемые ограничения и контрольные тесты;
- [`docs/GUI.md`](docs/GUI.md) — подробная история и концепция интерфейса;
- [`CLAUDE.md`](CLAUDE.md) — техническая карта живых/legacy-модулей и соглашения
  репозитория.

Идентификаторы в коде пишутся на английском; комментарии, docstring, логи и
сообщения об ошибках преимущественно на русском. Это соглашение проекта.
