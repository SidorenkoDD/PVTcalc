"""
Глобальные константы и допуски сходимости, используемые по всему движку.

TOL_TWO_PHASE_STABILITY_CONVERGENCE               — критерий сходимости пробных циклов TwoPhaseStabilityTest (см. _check_convergence).
TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION — критерий "тривиального решения" (K->1) в TwoPhaseStabilityTest.
TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE          — допуск бисекции Rachford-Rice в PhaseEquilibriumNewton.find_solve_bisection_v4.
TOL_TWO_PHASE_FLASH_CONVERGENCE                    — критерий сходимости по фугитивностям в PhaseEquilibriumNewton.check_convergence_ri.
TOL_TWO_PHASE_FLASH_TRIVIAL_SOLUTION               — критерий тривиального решения в PhaseEquilibriumNewton.check_trivial_solution.
TOL_SAT_PRESSURE                                   — допуск сходимости поиска давления насыщения (PhaseEnvelope/*).
CONSTANT_R                                         — универсальная газовая постоянная. Известный баг (CLAUDE.md Known Issues):
    округлено до 8.31, точное значение 8.314462618... — участвует во всех расчётах
    Z-фактора, молярного объёма и плотности по всему проекту. Не исправлено намеренно.
"""

TOL_TWO_PHASE_STABILITY_CONVERGENCE = 1e-12
TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION = 1e-4
TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE = 1e-13
TOL_TWO_PHASE_FLASH_CONVERGENCE = 1e-14
TOL_TWO_PHASE_FLASH_TRIVIAL_SOLUTION = 1e-09
TOL_SAT_PRESSURE = 1e-8
CONSTANT_R = 8.31