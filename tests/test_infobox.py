def test_infobox_classes():
    """Test `InfoBox` classes."""
    from aiidalab_qe.common.infobox import InfoBox

    custom_classes = ["custom-1", "custom-2 custom-3"]
    infobox = InfoBox(classes=custom_classes)
    assert all(
        css_class in infobox._dom_classes
        for css_class in (
            "info-box",
            "custom-1",
            "custom-2",
            "custom-3",
        )
    )


def test_in_app_guide():
    """Test `InAppGuide` class."""
    import ipywidgets as ipw

    from aiidalab_qe.common.guide_manager import guide_manager
    from aiidalab_qe.common.infobox import InAppGuide

    in_app_guide = InAppGuide(children=[ipw.HTML("Hello, World!")])
    assert all(
        css_class in in_app_guide._dom_classes
        for css_class in (
            "info-box",
            "in-app-guide",
        )
    )
    assert in_app_guide.children[0].value == "Hello, World!"

    assert in_app_guide.layout.display == "none"
    guide_manager.active_guide = (
        "Relaxation and electronic structure/Basic tutorial (relaxation, bands, PDOS)"
    )
    assert in_app_guide.layout.display == "flex"
    guide_manager.active_guide = "No guides"
    assert in_app_guide.layout.display == "none"

    guide_manager.active_guide = (
        "Relaxation and electronic structure/Basic tutorial (relaxation, bands, PDOS)"
    )
    in_app_guide = InAppGuide(identifier="guide-header")
    assert "You've activated the basic in-app guide" in in_app_guide.children[0].value


def test_plugin_guides_load_on_demand(monkeypatch):
    import aiidalab_qe.common.guide_manager as guide_manager_module

    calls = []
    monkeypatch.setattr(
        guide_manager_module,
        "get_entry_items",
        lambda *args: calls.append(args) or {},
    )

    manager = guide_manager_module.GuideManager()
    assert calls == []

    manager.load_plugin_guides()
    manager.load_plugin_guides()

    assert calls == [("aiidalab_qe.properties", "guides")]
