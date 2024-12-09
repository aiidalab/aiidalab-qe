from aiidalab_qe.common.infobox import InAppGuide, InfoBox


def test_infobox_classes():
    """Test `InfoBox` classes."""
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
    guide_class = "some_guide"
    in_app_guide = InAppGuide(guide_class=guide_class)
    assert all(
        css_class in in_app_guide._dom_classes
        for css_class in (
            "info-box",
            "in-app-guide",
            f"{guide_class}-guide-identifier",
        )
    )
