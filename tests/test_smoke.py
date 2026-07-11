"""Проверка, что pytest резолвит импорты calc_core.* без плясок с рабочей директорией."""


def test_calc_core_is_importable():
    from calc_core.Composition.Composition import Composition
    from calc_core.EOS.BaseEOS import EOSType

    assert Composition is not None
    assert EOSType.BRSEOS == "BRSEOS"
