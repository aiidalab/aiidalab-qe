from aiidalab_qe.common.infobox import FirstVisitBox, InAppGuide, InfoBox


def test_infobox_classes():
    """Test `InfoBox` class."""
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
    guide_id = "some_guide"
    in_app_guide = InAppGuide(guide_id=guide_id)
    assert all(
        css_class in in_app_guide._dom_classes
        for css_class in (
            "info-box",
            "in-app-guide",
            guide_id,
        )
    )


def test_first_visit_box():
    """Test `FirstVisitBox` class."""
    message = "This is a test"
    custom_classes = ["custom-1"]
    first_visit_box = FirstVisitBox(message=message, classes=custom_classes)
    assert all(
        css_class in first_visit_box._dom_classes
        for css_class in (
            "info-box",
            "first-visit-box",
            "custom-1",
        )
    )
    assert first_visit_box.children[0] == first_visit_box.message_box
    first_visit_box.close_button.click()
    assert first_visit_box.children[0] == first_visit_box.closing_message
    first_visit_box.undo_button.click()
    assert first_visit_box.children[0] == first_visit_box.message_box
    assert first_visit_box.message_box.children[0].value == message  # type: ignore
