# PVTcalc

## Project Overview

Самописный движок компонентного PVT-моделирования флюидов (нефть/газ/конденсат) на Python. Реализует полный расчётный цикл: описание состава флюида → расчёт равновесных свойств компонентов (в т.ч. псевдокомпонентов C7+) → уравнение состояния (EOS) → тест термодинамической стабильности → двухфазный флэш → стандартные PVT-эксперименты (CCE, DLE, сепарация, CVD) → построение фазовой диаграммы (давление насыщения, давление конденсации, критическая точка, фазовая огибающая).

Основное и единственное живое EOS — кубическое уравнение состояния **Брусиловского** (распространено в русской PVT-практике, обобщение Пенга-Робинсона/Соаве-Редлиха-Квонга). Legacy-реализации PREOS/SRKEOS удалены при чистке 2026-07-11 — см. [Live vs Legacy Modules](#live-vs-legacy-modules-важно).

Автор — практикующий инженер-разработчик модели, ведёт разработку в одиночку (иногда с правками коллеги "Дениса" — см. файлы `*_den_v.py`). По словам автора, текущая версия расчётов сходится с результатами коммерческого симулятора **PVTSim**. Явного регрессионного теста, фиксирующего это сравнение в репозитории, пока нет (см. [Known Issues](#known-issues--bugs) и рекомендации).

**Язык кода**: идентификаторы (классы/функции/переменные) — на английском, комментарии/docstring/сообщения об ошибках/логи — **на русском**. Это осознанная и последовательная конвенция автора, не ошибка — не "исправлять" при рефакторинге.

## Environment & Running

- Python **3.11+** обязателен (используется `enum.StrEnum` в `calc_core/EOS/BaseEOS.py`). Фактически код гонялся под 3.12 и 3.13 (видно по `__pycache__/*.cpython-312.pyc` и `*.cpython-313.pyc`).
- **Зависимости и сборка — `pyproject.toml`** (с 2026-07-11): `[project.dependencies]` — рантайм (`numpy`, `pandas`, `scipy`, `matplotlib`, `joblib`, `openpyxl`), `[project.optional-dependencies]` — `dev` (pytest) и `notebook` (`ipykernel`, `jupyterlab`) отдельно. `[build-system]` (setuptools) настроен и рабочий — `pip install -e ".[dev]"` устанавливает пакет в editable-режиме: правки файлов в `calc_core/` подхватываются сразу, без переустановки.
  - Distribution-имя пакета (то, что видит `pip install`) — `pvtcalc`; импортируемое имя внутри Python — `calc_core` (`from calc_core.VLE.Flash import Flash`). Изначально (тестовая сборка 2026-07-11) импортируемым именем был `_src` — переименовано в `calc_core` в тот же день по решению автора (совпадает с именем корневой папки проекта на его стороне). Проверено тестовым `.whl`-билдом в изолированном venv до и после переименования — оба раза работает без каких-либо repo-relative путей.
  - До этой правки пакета не существовало вообще (`pip install -e .` было невозможно, все внутренние импорты держались на том, что `test_notebook.ipynb` физически лежит в корне репозитория и Jupyter добавляет cwd в `sys.path`). После `pip install -e ".[dev]"` эта случайность больше не нужна — импорты `calc_core.*` резолвятся из любой рабочей директории.
- Все подпакеты `calc_core/**` теперь имеют `__init__.py` (раньше 6 из 9 были неявными namespace-пакетами без него — работало благодаря PEP 420, но мешало однозначному discovery для сборки). Файлы пустые, кроме `calc_core/PlusComponents/__init__.py`, где по ошибке лежит код верхнего уровня, а не re-export. Никакого `__all__`/переэкспорта нет — всё импортируется напрямую из конкретных модулей.

## Directory Structure

```
PVTcalc/
├── CLAUDE.md
├── models.json                 # сохранённые "снэпшоты" моделей флюидов (состав + свойства + результаты расчётов)
├── test_notebook.ipynb         # ЕДИНСТВЕННЫЙ актуальный интеграционный сценарий / рабочий скретчпад автора
├── gui/                        # DearPyGui UI: state/services/view + session/job orchestration
├── tests/                      # unit/regression/headless GUI tests
└── calc_core/
    ├── Composition/            # представление состава флюида (Composition.py)
    ├── CompositionalModel/     # верхнеуровневый фасад Composition+EOS+Flash (CompositionalModel.py)
    ├── EOS/                    # уравнения состояния: Брусиловский (текущий), PR/SRK (legacy, удалены)
    ├── Experiments/            # стандартные PVT-эксперименты: CCE, DLE, сепарация, CVD
    ├── PhaseEnvelope/          # давление насыщения/конденсации, критическая точка, фазовая огибающая
    ├── PhaseStability/         # тест стабильности Михельсена, идентификатор фазы
    ├── PlusComponents/         # корреляции свойств псевдокомпонентов C7+ (по одному файлу на свойство)
    ├── Utils/                  # константы, БД компонентов (DB.json), EOS-парам. БД, вязкость, persist-слой
    └── VLE/                    # решатель двухфазного равновесия (Flash, Ньютон по фугитивностям)
```

## Live vs Legacy Modules (важно!)

Репозиторий проходил через незавершённый рефакторинг/переименование пакета (было что-то вроде `calculations.*`, стало `_src.*`, затем `_src.*` переименован в `calc_core.*` — см. ниже). Часть файлов не была ни обновлена, ни удалена — они **импортировали несуществующий пакет `calculations.*` и не запускались в принципе** (`ImportError` при первом же импорте). **2026-07-11**: 24 таких файла удалены из репозитория (см. историю коммитов) после трассировки импортов по всему проекту, включая ноутбук — ни один не имел живых потребителей. (`calc_core/Utils/BaseClasses.py` изначально тоже попал под удаление, но был восстановлен в тот же день: его `calculations.*`-импорт был закомментирован, а сам файл реально используется — `JsonDBReader.py` наследует от него `Reader`.)

**2026-07-11 (тем же днём)**: после удаления старых версий суффиксы `V2`/`V3`/`2` в именах живых файлов потеряли смысл (различать стало не с чем) и были убраны через `git mv` — с обновлением всех импортов в `_src/` (тогдашнее имя пакета) и в `test_notebook.ipynb`. Старое имя `DB.json` (471 строка, legacy-схема) было удалено, а `DB_V2.json` переименован в `DB.json` (961 строка, актуальная схема с `available_components`/`sequence_number`/`carbon_flag`). Класс-имена (`Composition`, `BrusilovskiyEOS`, `TwoPhaseStabilityTest` и т.д.) суффиксов никогда не содержали — менялись только имена файлов и путей импорта. **Если встречаешь в старых заметках/памяти упоминания `CompositionV2.py`, `BrusilovskiyEOSV2.py`, `TwoPhaseStabilityTestV3.py`, `DLE2.py`, `CCE2.py`, `SeparatorTest2.py`, `Export2.py`, `Results2.py`, `PhaseEquilibriumNewtonV2.py`, `FluidPropertiesCalculatorV2.py`, `BRS_EOS_DB_V2.py`, `CompositionalModelV2.py`, `DB_V2.json` — это устаревшие имена, актуальные файлы называются так же, но без суффикса.**

**2026-07-11 (позже тем же днём)**: пакет `_src` целиком переименован в `calc_core` (`git mv _src calc_core` + замена `_src.` → `calc_core.` во всех импортах — `calc_core/**`, `tests/`, `test_notebook.ipynb`, `pyproject.toml`, `docs/diagrams/module-interactions.html`). Мотивация — совпадение с именем корневой папки проекта у автора; до этого пакет назывался `_src` начиная с самой первой чистки того же дня. **Если встречаешь в старых заметках/памяти упоминания `_src.*` или путей `_src/...` — это устаревшее имя, актуальное имя пакета — `calc_core`.**

Ниже — актуальная таблица "что живое" по каждой области. **При любой задаче открывай и правь только левую колонку.**

| Область | Канонический (живой) модуль | Не использовать (legacy / сломано) |
|---|---|---|
| Состав флюида | `calc_core/Composition/Composition.py` (`Composition`) | — (старые `Composition.py`/`component.py`/`bips.py` первого поколения удалены) |
| БД компонентов | `calc_core/Utils/DB.json` (через `JsonDBReader`) | — (старая схема `DB.json` удалена, на её место переименован бывший `DB_V2.json`) |
| EOS | `calc_core/EOS/BrusilovskiyEOS.py` (`BrusilovskiyEOS`) | `calc_core/EOS/RootChooser.py` — **сохранён намеренно по решению автора**, хотя импортирует `calculations.*` и не запускается; `MixingRule.py` — не сломан, но не используется (см. категорию "осиротевший код" ниже) |
| Тест стабильности | `calc_core/PhaseStability/TwoPhaseStabilityTest.py` (`TwoPhaseStabilityTest`) | — (старые версии/`BasePhaseStability.py` удалены) |
| Двухфазный флэш | `calc_core/VLE/Flash.py` + `calc_core/VLE/PhaseEquilibriumNewton.py` | — (`PhaseEquilibrium.py` удалён) |
| Модель-фасад | `calc_core/CompositionalModel/CompositionalModel.py` | — (`Variant.py` удалён) |
| PVT-эксперименты | `calc_core/Experiments/DLE.py`, `CCE.py`, `SeparatorTest.py` | — (`StandardSeparation.py`/`ExperimentsFacade.py` удалены) |
| Persist-слой | `calc_core/Utils/Export.py` (`ModelJSONDB`), `calc_core/Utils/Results.py` (`ResultStore`) | — (`ResultsViewer.py` удалён) |
| Фазовая огибающая | `calc_core/PhaseEnvelope/PhaseEnvelopeSuccessiveSubstitution.py` (`PhaseEnvelopeSSM`, основная реализация) + `PhaseEnvelopeNewton.py` (`PhaseEnvelopeNewton`, точнее в отдельной точке, но менее устойчив у критической точки) + `BubblePointPressure.py`/`DewPressure.py` (`BubblePointCalculator`/`DewPointCalculator`, Ньютон по производным летучести) + `PhaseEnvelopeGrid.py` (`PhaseEnvelopeGrid`, сеточная затравка нижней ветки SSM) + `CriticalPoint.py` (`CriticalPointCalculator`, минимизация зазора bubble/dew) | `PhaseDiagram_v4.py`, `PhaseEnvelope.py` — удалены 2026-07-15 при выборе канонической реализации (см. ниже) |

Единственный оставшийся файл с импортом несуществующего `calculations.*` — **`calc_core/EOS/RootChooser.py`**, сохранён по явному решению автора (не удалять и не трогать без отдельного запроса).

**2026-07-15**: папка `calc_core/PhaseDiagram/` переименована в `calc_core/PhaseEnvelope/` (`git mv`), одновременно выбрана каноническая реализация фазовой огибающей/давления насыщения/критической точки из ранее конкурировавших вариантов (см. таблицу выше и обновлённый раздел [Зона активной разработки](#зона-активной-ещё-не-устоявшейся-разработки)). Удалены `PhaseDiagram_v4.py` (legacy, импортировал несуществующий `calculations.*`) и `PhaseEnvelope.py` (старый wrapper, не имел живых потребителей в `calc_core/`); `PhaseEnvelopeFromStability.py` переименован в `PhaseEnvelopeGrid.py` (класс `PhaseEnvelopeFromStability` → `PhaseEnvelopeGrid`). При этом в папке намеренно остались НЕ канонические, но живые файлы: `new_methodv2.py::SaturationPressure` (используется в `DLE.py`/`CCE.py`/`SeparatorTest.py` — авто-определяет bubble/dew через тест стабильности, замена на `BubblePointCalculator` физически не эквивалентна для газоконденсатных систем, решение отложено до связки с сервисным слоем/`CompositionalModel`) и `SaturationPressure_den_v.py`/`CriticalProperties_den_v.py` (нужны только `CVD_den_v.py`, см. [осиротевший код](#осиротевшийнезавершённый-код-не-путать-с-legacy--технически-рабочий-но-нужно-решение-автора) — путь у обоих файлов при переименовании папки тоже поменялся). **Если встречаешь в старых заметках/памяти упоминания `calc_core/PhaseDiagram/...` или `PhaseEnvelopeFromStability` — это устаревшие пути/имя, актуальная папка — `calc_core/PhaseEnvelope/`, актуальный класс — `PhaseEnvelopeGrid`.**

### Осиротевший/незавершённый код (не путать с legacy — технически рабочий, но нужно решение автора)

- `calc_core/Experiments/CVD_den_v.py` — **сейчас сломан**, но не из-за `calculations.*`: импортирует `SaturationPointCalculator` из `SaturationPressure_den_v.py`, а там определён только класс `DewPointCalculator` — такого имени в файле нет → `ImportError`. Нигде не импортируется (в т.ч. не используется в ноутбуке, несмотря на то что CVD как эксперимент нигде больше не реализован). Похоже на недописанную интеграцию, а не на легаси — вероятно, стоит доделать (поправить импорт), а не удалять. Автор подтвердил: этот модуль понадобится в будущем, пока не трогаем.
- `calc_core/PhaseEnvelope/CriticalProperties_den_v.py` — рабочий, отдельный алгоритм критической точки (критерий касательной плоскости), но **единственный потребитель — сломанный `CVD_den_v.py`**; в ноутбуке не используется (там используется `CriticalPoint.py::CriticalPointCalculator`). Сейчас фактически недостижим, но исходно рабочий код "Дениса" — решение об удалении/интеграции за автором.
- `calc_core/Experiments/BaseExperiment.py` — импортируется только сломанным `CVD_den_v.py`; судьба привязана к решению по нему.
- `calc_core/EOS/MixingRule.py` (`MixingRule`, `ClassicMixingRule`, `HuronVidalMixingRule`, `MixingRuleFactory`) — самодостаточный, не сломан, но нигде не используется — похоже на заготовку под будущую поддержку разных правил смешения (сейчас смешение зашито прямо в `BrusilovskiyEOS.py`).
- `calc_core/Utils/Import.py` (`DBModelImport`) — класс с пустым `__init__`, нигде не вызывается. Автор подтвердил: этот модуль точно понадобится — по смыслу должен стать зеркалом `Export.py::ModelJSONDB` (тот пишет `Composition` → `models.json`, этот должен читать `models.json` → `Composition`/список моделей).
- `calc_core/PhaseStability/PhaseStabilityPhasade.py` — пустой файл (0 байт), нигде не используется.
- `calc_core/Composition/component_tests.py`, `calc_core/Composition/composition_v2_test.py` — ручные smoke-скрипты с относительными импортами (`from component import Component`), работают только если `calc_core/Composition/` — рабочая директория; не pytest, не CI. Не легаси в смысле "сломано", но и не часть рабочего пайплайна.

Решение по всей этой категории отложено намеренно — держать в памяти как "недоработанное", не как "мёртвое", до отдельного захода.

### Зона активной, ещё не устоявшейся разработки

**2026-07-15**: неопределённость с канонической реализацией `calc_core/PhaseEnvelope/` (бывшая `PhaseDiagram/`) разрешена — см. строку "Фазовая огибающая" в таблице "что живое" выше. Кратко: основной путь — `PhaseEnvelopeSSM` (устойчив у критической точки), точечно — `PhaseEnvelopeNewton`/`BubblePointCalculator`/`DewPointCalculator` (точнее, но менее устойчивы вблизи критической точки), затравка нижней ветки SSM — `PhaseEnvelopeGrid`, критическая точка — `CriticalPoint.py::CriticalPointCalculator`. `new_methodv2.py`, `SaturationPressure_den_v.py`, `CriticalProperties_den_v.py` в папке остались, но канонической реализацией не являются (живые зависимости/отложенная работа, см. заметку 2026-07-15 выше и раздел "осиротевший код").

Остаётся нерешённым (не входило в эту задачу): в `test_notebook.ipynb` всё ещё существуют **только внутри ячеек**, не вынесенные в `calc_core/`, ещё как минимум две реализации критической точки/фазовой диаграммы — `CriticalPointCalculator`/`CriticalPointCalculatorParrallel` (сеточный поиск по метрике близости S→1) и `PhaseDiagramCalculator` в стиле VBA-бисекции (секция "#### From VBA"). Если ставится задача, связанная именно с этими двумя — уточнить у пользователя, а не угадывать.

## Calculation Flow

> Визуальная версия того, что описано ниже текстом (плюс слоистая карта пакетов `calc_core/` и карта зоны дублирования `PhaseDiagram/`) — см. [`docs/diagrams/module-interactions.html`](docs/diagrams/module-interactions.html) (откройте в браузере).

Основной живой путь расчёта одного флэша:

```
Composition(zi, T_res, eos_name)                    calc_core/Composition/Composition.py
  └─ evaluate_composition_data(...)                  Tc/Pc/ω/Vc/shift/Kw на компонент
       ├─ PlusComponentProperties (для C7+)           calc_core/PlusComponents/*
       └─ BRS_EOS_DB (параметры EOS + BIP)            calc_core/Utils/BRS_EOS_DB.py

Conditions(p_bar, t_C)                                calc_core/Utils/Conditions.py  (t_K = t_C + 273.15)

Flash(composition, conditions).calculate()            calc_core/VLE/Flash.py
  ├─ TwoPhaseStabilityTest(...)                        calc_core/PhaseStability/TwoPhaseStabilityTest.py
  │     └─ BrusilovskiyEOS(...) — фугитивности          calc_core/EOS/BrusilovskiyEOS.py
  ├─ если 2 фазы: PhaseEquilibriumNewton(...)           calc_core/VLE/PhaseEquilibriumNewton.py
  │     Rachford-Rice (Ньютон + бисекция) → уточнение K по фугитивностям (аналитический якобиан)
  ├─ FluidPropertiesCalculator(...) на каждую фазу       calc_core/Utils/FluidPropertiesCalculator.py
  │     молярная масса, молярный объём, плотность, Z со сдвигом, вязкость (LBC, calc_core/Utils/Viscosity.py)
  └─ → FlashResult(vapor: PhaseState, liquid: PhaseState, is_two_phase)   calc_core/VLE/FlashResult.py
```

Однофазная ветка `Flash.calculate()` сохраняет историческую раскладку долей фаз (весь состав в liquid), необходимую для совместимости DLE, но с 2026-07-18 отдельно возвращает `FlashResult.phase_type`: `vapor`, `liquid` либо `ambiguous`. Тип определяет `calc_core/PhaseStability/PhaseIdentificator.py`; старое ошибочное имя класса оставлено как alias.

Эксперименты — надстройка над `Flash` + `SaturationPressure`, работают по ступеням (давления или пар давление/температура), передавая состав жидкости с предыдущей ступени на следующую через `composition.new_composition(...)`:

```
DLE.DLE.calculate() / SeparatorTest.SeparatorTest.calculate() / CCE.CCE.calculate()
  → new_methodv2.SaturationPressure(...) — находим точку насыщения
  → цикл по ступеням: Flash(...) на составе с предыдущей ступени
  → Bo, Rs и т.п. из серии FlashResult
```

## Data Model & Databases

- Состав флюида — **не dataclass**, а "колоночная" структура словарей: `composition: Dict[str, float]` (мольные доли) + `composition_data: Dict[property_name, Dict[component_name, value]]`. Объекта отдельного компонента в живом коде нет (старый `component.py::Component` был legacy и удалён при чистке 2026-07-11).
- Статические свойства чистых компонентов (молярная масса, Tb, критические параметры, ацентрический фактор, shift-параметр, BIP, Kw) — в `calc_core/Utils/DB.json`, читаются через `calc_core/Utils/JsonDBReader.py` (наследует абстрактный `Reader` из `calc_core/Utils/BaseClasses.py` — единственный живой класс в этом файле, остальные его абстракции (`PhaseStabilityTest`, `Calculator`, `CalculationModule`) реальным кодом не используются; сам `JsonDBReader` умеет искать файл по нескольким кандидатным путям — cwd, родительские директории, `$DB_PATH`, каталог `sys.argv[0]`; это позволяет запускать код из разных рабочих директорий, но одновременно скрывает ошибку, если файл не найден там, где ожидалось).
- Свойства псевдокомпонентов C7+ считаются на лету корреляциями в `calc_core/PlusComponents/` (по одному файлу на свойство: `CriticalTemperature.py`, `CriticalPressure.py`, `CriticalVolume.py`, `AcentricFactor.py`, `ShiftParameter.py`, `kWatson.py`), диспетчеризуются по строковым ключам метода (например `{'critical_temperature': 'pedersen', 'critical_pressure': 'rizari_daubert', ...}`). Это и есть результат рефакторинга "компоненты разнесены по модулям и в отдельной папке" — единственный источник этой логики.
- `models.json` (корень репозитория) — персистентный "снэпшот"-стор: полный состав + рассчитанные свойства + история результатов флэша для именованных моделей (сейчас — `KRSNL_PVTSIM`/"KRASNOLEN_PVTSIM_COMPARISON" и `PRRZLM_MDT_TEST`). Пишется через `calc_core/Utils/Export.py::ModelJSONDB`, читается через `Composition.from_db(path)`. **Важно**: это не независимый эталонный датасет PVTSim — сюда пишутся результаты расчётов самого движка PVTcalc, а название `KRSNL_PVTSIM` означает лишь то, что этот конкретный состав когда-то сверялся с PVTSim во внешнем инструменте, но само численное сравнение в репозитории не зафиксировано.

## Conventions

- Комментарии/docstring/сообщения об ошибках — на русском, идентификаторы — на английском. Исключение — `shift_parametr` (транслитерация "parameter"), устойчиво используется в `BrusilovskiyEOS.py`, `FluidPropertiesCalculator.py` — не переименовывать без явного запроса пользователя, это breaking change по многим файлам.
- **Суффиксы версий в именах файлов больше не используются** (сняты 2026-07-11, см. [Live vs Legacy Modules](#live-vs-legacy-modules-важно)). При добавлении новой реализации предпочитать **замену** старой (или явную договорённость с пользователем о новой версии/суффиксе), а не создание параллельной ветки без необходимости — в проекте и так есть зона с 3+ реализациями одной задачи (см. [Зона активной разработки](#зона-активной-ещё-не-устоявшейся-разработки)). Папка `PhaseDiagram/` — единственное намеренное исключение, там суффиксы (`new_methodv2.py` и т.п.) сохранены до отдельной работы над этой зоной.
- `test_notebook.ipynb` — не второстепенный файл, а фактический единственный интеграционный тест и источник истины о том, какой пайплайн реально работает. При изменении любого модуля в `calc_core/` полезно свериться с соответствующей секцией ноутбука (секции подписаны markdown-заголовками на русском: "Расчет флеша", "Давление насыщения...", "CCE в параллель" и т.д.), что он не сломан логически (формального прогона/CI нет).
- Единицы измерения: давление в барах, температура на входе в `Conditions(p, t)` подаётся в °C и внутри конвертируется в Кельвины. При работе с фазовыми огибающими (`calc_core/PhaseEnvelope/`) явно проверять, какие единицы ожидает конкретный модуль (исторически там встречались расхождения — legacy-модуль на МПа удалён 2026-07-15, но привычку проверять стоит сохранить).
- Численные изменения ядра должны подтверждаться регрессионными тестами. Исправления 2026-07-18 (`273.15`, точная `R`) выполнены по прямому запросу автора убрать существующие ошибки; baseline Flash обновлён по фактическим новым значениям.

## Logging

Единый подход к логированию (введён 2026-07-11, см. `calc_core/Utils/Logging.py`):

- Библиотечный код (`calc_core/**`) **никогда** не вызывает `logging.basicConfig()` и не навешивает handlers сам — только объявляет `logger = logging.getLogger(__name__)` в начале модуля и пишет через него. `__name__` внутри пакета даёт естественную иерархию (`calc_core.VLE.Flash`, `calc_core.EOS.BrusilovskiyEOS` и т.д.), корень которой — логгер `calc_core`.
- Настройка вывода (куда писать, какой уровень) — дело того, кто использует движок: вызвать `calc_core.Utils.Logging.configure_logging(level=logging.INFO)` один раз в начале ноутбука/скрипта/теста. Функция идемпотентна.
- До этой правки три модуля в `PhaseDiagram/` сами вызывали `logging.basicConfig()`/навешивали handler как побочный эффект импорта (конфликтовали друг с другом, порядок импорта решал итоговую конфигурацию) — убрано. `CriticalPoint.py::CriticalPointCalculator` сохранил параметр `verbose` в сигнатуре для обратной совместимости, но он больше ничего не настраивает.
- Основной расчётный путь (`Flash`, `PhaseEquilibriumNewton`, `TwoPhaseStabilityTest`, `BrusilovskiyEOS`, `CompositionalModel`, `DLE`/`CCE`/`SeparatorTest`) до этой правки не логировал вообще ничего — теперь логирует на уровнях `DEBUG` (детали итераций/выбор корня EOS) и `INFO` (стабильно/нестабильно, P_sat, прогресс по ступеням). Полезно при отладке нестабильностей около критической точки — например, `BrusilovskiyEOS.calc_eos()` на DEBUG явно показывает, когда кубическое уравнение даёт несколько практически совпадающих действительных корней.

## Known Issues / Bugs

Исправлено 2026-07-18: смещение температуры `273.15`, стандартная температура `293.15 K`, точная `CONSTANT_R`, режим `wb` у `ResultStore.save()`, канонический `LengthMismatchError` (старое имя — совместимый alias), возвращаемый тип фазы и каноническое имя `PhaseIdentificator`, правила `.gitignore` для `*.pyc`/`.DS_Store`.

Остаются намеренно вне живого пути:

- `calc_core/EOS/RootChooser.py` всё ещё импортирует удалённый `calculations.*`; файл сохранён по решению автора и не используется.
- `calc_core/Experiments/CVD_den_v.py` ожидает старый API `SaturationPointCalculator`; интеграция CVD отложена до отдельной предметной задачи.
- Отмена GUI-задачи не прерывает итерации внутри солвера, а только отбрасывает результат; для настоящей отмены нужен cancellation token в расчётном API.
- `gui_settings.json` пока не подключён к константам движка; до появления явного config object это только сохранённый профиль UI.

Список намеченной (пока не сделанной) работы по этим и другим пунктам — архитектура, кодовая база, тестирование, защита/валидация — ведётся отдельно в [`docs/BACKLOG.md`](docs/BACKLOG.md), с приоритетами P0/P1/P2. Этот раздел фиксирует факты о текущем состоянии кода, `BACKLOG.md` — план работы по ним.

## Recent Development Focus

По git log на момент написания (ветка `correct_calc`): последние коммиты сосредоточены на (1) идентификации фазы в однофазном состоянии выше крикондентермы и вблизи критической точки — область ещё нестабильна, и (2) разнесении корреляций свойств псевдокомпонентов C7+ по отдельным модулям в `calc_core/PlusComponents/`. 2026-07-11 добавлена ветка `cleanup/remove-dead-legacy-modules`: удалена мёртвая `calculations.*`-легаси и сняты избыточные суффиксы версий (`V2`/`V3`/`2`) в именах живых файлов (кроме папки `PhaseDiagram/`, оставленной под отдельную работу). Это следует держать в уме как текущий фокус активной разработки при планировании новых задач.

2026-07-15: заведён `docs/BACKLOG.md` (приоритизированный P0/P1/P2 план развития — архитектура/кодовая база/тестирование/защита, см. ссылку в Known Issues) и выполнен первый P0-пункт из него: выбрана каноническая реализация фазовой огибающей, `calc_core/PhaseDiagram/` переименована в `calc_core/PhaseEnvelope/` с удалением конкурирующих legacy-реализаций (см. таблицу "что живое" и раздел "Зона активной разработки" выше).

2026-07-18: GUI прошёл стабилизационный рефакторинг без добавления предметных функций. Персистентность атомарна, редактируемые модели имеют dirty/save lifecycle, фоновые задачи адресуются через стабильный `NodeRef`, сессия v3 вынесена в `workspace_codec`, orchestration — в `calculation_coordinator`; добавлены headless DPG и persistence-тесты. Полное описание — `docs/GUI.md`, СТАТУС-6.

2026-07-18 (следующий блок): монолитный `gui/view/app.py` декомпозирован на composition root (около 500 строк) и семь предметных `gui/view/*_view.py`: Projects, Composition, Flash/Compare, Experiments, Envelope, общие dialogs и Workspace tree/tabs. Поведение не менялось; headless smoke теперь собирает вкладки всех основных типов. Детали — `docs/GUI.md`, СТАТУС-7.
