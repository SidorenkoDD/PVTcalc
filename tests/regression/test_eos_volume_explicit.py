"""Regression-тест на (T,V,n)-нативную оценку BrusilovskiyEOS.eval_pv_lnphi/eval_at_volume.

Новая явная (без решения кубического уравнения) формула P(T,V,n)/ln(phi_i)(T,V,n)
добавлена для метода Михельсена определения критической точки
(`calc_core/PhaseDiagram/CriticalPointMichelsen.py`, работает в (T,V), а не в
(T,P) — см. CLAUDE.md). Формула выведена алгебраической подстановкой Z=PV/(RT)
в P=RT/(V-b)-a/((V+c)(V+d)) и сверена почленно с E0/E1/E2 из `_calc_roots_eos`
(см. докстринг `eval_pv_lnphi`) — этот тест проверяет то же самое численно,
на реальных составах: берём Z, найденный существующим (T,P)-путём (`calc_eos()`),
считаем V=ZRT/P и проверяем, что `eval_at_volume(V)` воспроизводит то же P и
тот же вектор ln(phi_i).
"""

import numpy as np
import pytest

from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.Utils.Constants import CONSTANT_R

REL = 1e-8


def _check_volume_explicit_matches_pt_path(composition, p: float, t: float):
    composition.T = t
    eos = BrusilovskiyEOS(composition, p, t)
    z_chosen, _ = eos.calc_eos()

    # V во ВНУТРЕННИХ единицах EOS (те же, что у composition_data['a'|'b'|'c'|'d']),
    # не см3/моль FluidPropertiesCalculator — см. докстринг eval_pv_lnphi.
    v = CONSTANT_R * t * z_chosen / p

    p_from_v, ln_phi_from_v = eos.eval_at_volume(v)
    ln_phi_from_pt = eos.get_fugacity_coef_vector_by_root(z_chosen)

    assert p_from_v == pytest.approx(p, rel=REL)
    np.testing.assert_allclose(ln_phi_from_v, ln_phi_from_pt, rtol=REL, atol=1e-10)


def test_volume_explicit_matches_pt_path_krsnl_single_phase(krsnl_composition):
    """KRSNL_PVTSIM, P=250 бар, T=390°C — та же точка, что test_flash.py::test_krsnl_single_phase_above_p_sat."""
    _check_volume_explicit_matches_pt_path(krsnl_composition, p=250.0, t=390.0 + 273.15)


def test_volume_explicit_matches_pt_path_krsnl_two_phase_region(krsnl_composition):
    """KRSNL_PVTSIM, P=100 бар, T=110°C — та же точка, что test_flash.py::test_krsnl_two_phase_below_p_sat
    (двухфазная область; здесь проверяется EOS смеси как таковой, без сборки Flash)."""
    _check_volume_explicit_matches_pt_path(krsnl_composition, p=100.0, t=110.0 + 273.15)


def test_volume_explicit_matches_pt_path_przlm(przlm_composition):
    """PRRZLM_MDT_TEST, P=200 бар, T=100°C — та же точка, что test_flash.py::test_przlm_single_phase."""
    _check_volume_explicit_matches_pt_path(przlm_composition, p=200.0, t=100.0 + 273.15)
