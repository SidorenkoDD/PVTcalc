"""Проверка, что pytest резолвит импорты _src.* без плясок с рабочей директорией."""


def test_src_is_importable():
    from _src.Composition.Composition import Composition
    from _src.EOS.BaseEOS import EOSType

    assert Composition is not None
    assert EOSType.BRSEOS == "BRSEOS"
