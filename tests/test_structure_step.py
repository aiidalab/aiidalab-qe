import builtins

import aiidalab_widgets_base as awb
from aiidalab_qe.app.structure.step import _get_optional_cdxml_importer


def test_optional_cdxml_importer_is_used_when_available(monkeypatch):
    class DummyCdxmlUploadWidget:
        def __init__(self, title):
            self.title = title

    monkeypatch.setattr(awb, "CdxmlUploadWidget", DummyCdxmlUploadWidget, raising=False)

    importer = _get_optional_cdxml_importer()

    assert isinstance(importer, DummyCdxmlUploadWidget)
    assert importer.title == "CDXML"


def test_optional_cdxml_importer_is_skipped_when_unavailable(monkeypatch):
    original_import = builtins.__import__

    def import_without_cdxml(name, globals_=None, locals_=None, fromlist=(), level=0):
        if name == "aiidalab_widgets_base" and "CdxmlUploadWidget" in fromlist:
            raise ImportError("CdxmlUploadWidget is unavailable")
        return original_import(name, globals_, locals_, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", import_without_cdxml)

    assert _get_optional_cdxml_importer() is None
