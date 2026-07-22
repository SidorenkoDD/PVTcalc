# Корреляции псевдокомпонентов C7+

Статус на 2026-07-21: все корреляции, зарегистрированные в
`calc_core/PlusComponents`, имеют внутренние контрольные значения и негативные
тесты. Это проверка реализации формул, единиц и диспетчеризации, но не замена
сверке с первичными публикациями для конкретного флюида.

## Публичный путь

`Composition.evaluate_composition_data()` создаёт `PlusComponentProperties` для
каждой компоненты с `c7_plus_flag=True`. Агрегатор считает свойства в порядке:

```text
Kw → critical_temperature → critical_pressure
   → acentric_factor → critical_volume → shift_parameter
```

Порядок нужен корреляциям, использующим уже вычисленные свойства. Явно заданные
в конструкторе `Tc`, `Pc`, `Kw` и `af` считаются входными данными и не
перезаписываются `calculate_all()`.

Базовые единицы:

- `M` — g/mol;
- `gamma` — безразмерная относительная плотность;
- `Tb`, `Tc` — K;
- `Pc` — bar;
- `Vc` — cm³/mol;
- `Kw`, `acentric_factor`, `shift_parameter` — безразмерные.

## Зарегистрированные методы

| Свойство | Выход | Требуемые входы → методы |
|---|---|---|
| Critical temperature | K | `gamma,Tb` → Roess, Nokey, Cavett, Kesler–Lee, Sim–Daubert, Riazi–Daubert; `gamma,M` → Pedersen, Standing, Mogoulas–Tassios; `Tb` → Twu; `gamma,Tb,M` → Watansiri–Owens–Starling |
| Critical pressure | bar | `gamma,Tb` → Kesler–Lee, Riazi–Daubert, Cavett, Sim–Daubert; `gamma,M` → Pedersen, Standing, Mogoulas–Tassios |
| Acentric factor | 1 | `Pc,Tc,Tb` → Edmister; `M` → Riazi–Al-Sahhaf; `Pc,Tc,Tb,Kw` → Kesler–Lee; `gamma,M` → Mogoulas–Tassios |
| Critical volume | cm³/mol | `gamma,M` → Hall–Yarborough, Lohrenz; `gamma,Tb` → Riazi–Daubert; `af,Pc,Tc` → Reid |
| Watson factor | 1 | `gamma,Tb` → K Watson; `gamma,M` → Riazi–Daubert |
| Shift parameter | 1 | `M,Kw` → Jhaveri–Youngren; `af` → SRK или PR; пустой метод → `0.0` |

Незавершённые Twu/Watansiri для `Pc`, Riedel для `Vc` и BRS shift не
зарегистрированы и через публичный агрегатор недоступны.

## Исполняемые ограничения

- `M`, `gamma`, `Tb`, а также заданные `Tc`, `Pc`, `Kw` должны быть конечными
  положительными числами; `af` — конечным числом.
- Рассчитанные `Tc`, `Pc`, `Vc`, `Kw` должны быть конечными и положительными.
- Standing: `M > 71.2 g/mol` для `Tc`, `M > 61.1 g/mol` для `Pc`.
- Jhaveri–Youngren выбирает коэффициенты по диапазонам `Kw`: aromatic
  `(8.5, 11)`, naphthenic `(11, 12.5)`, paraffinic `(12.5, 13.5)`; вне них и
  ровно на границах текущая реализация возвращает `0.0`.
- Неизвестное свойство/имя метода, отсутствующая зависимость или
  незарегистрированный TODO дают явную ошибку.

Точные опубликованные диапазоны применимости большинства исторических
корреляций в репозитории пока не сохранены. Поэтому интерфейс не должен
автоматически рекомендовать альтернативный метод только потому, что он
зарегистрирован: для этого нужен отдельный каталог источников и диапазонов.

## Контрольный псевдокомпонент

Матрица тестов использует `M=200 g/mol`, `gamma=0.85`, `Tb=600 K`; для методов
с зависимостями также `Kw=11.5`, `Tc=750 K`, `Pc=20 bar`, `af=0.6`.
Проверяются все 32 зарегистрированных метода. Отдельный сквозной тест закрепляет
используемую проектами конфигурацию Pedersen / Riazi–Daubert /
Riazi–Al-Sahhaf / Hall–Yarborough / K Watson / Jhaveri–Youngren.

В R1.5 исправлены:

- Reid `Vc`: добавлен коэффициент перевода `R` из J/(mol·K) при `Pc` в bar —
  прежний результат был ровно в десять раз меньше;
- Standing `Pc`: нижняя граница исправлена с `53.7` на `61.1 g/mol`, которую
  фактически требует `log10(M-61.1)`;
- `calculate_all()`: заранее переданные `Tc/Pc/Kw/af` больше не
  перезаписываются корреляциями.

Исполняемая спецификация: `tests/unit/test_plus_components.py`.
