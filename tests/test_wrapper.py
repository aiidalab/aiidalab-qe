from aiidalab_qe.app.wrapper import AppWrapperContoller, AppWrapperModel, AppWrapperView


class TestWrapper:
    def test_enable_controls(self):
        """Test enable_controls method."""
        self._instansiate_mvc_components()
        assert self.view.guide_toggle.disabled is True
        assert self.view.about_toggle.disabled is True
        self.controller.enable_controls()
        assert self.view.guide_toggle.disabled is False
        assert self.view.about_toggle.disabled is False

    def test_guide_toggle(self):
        """Test guide_toggle method."""
        self._instansiate_mvc_components()
        self.controller.enable_controls()
        self.controller._on_guide_toggle({"new": True})
        self._assert_guide_is_on()
        self.controller._on_guide_toggle({"new": False})
        self._assert_no_info()

    def test_about_toggle(self):
        """Test about_toggle method."""
        self._instansiate_mvc_components()
        self.controller.enable_controls()
        self.controller._on_about_toggle({"new": True})
        self._assert_about_is_on()
        self.controller._on_about_toggle({"new": False})
        self._assert_no_info()

    def test_toggle_switch(self):
        """Test toggle_switch method."""
        self._instansiate_mvc_components()
        self.controller.enable_controls()
        self._assert_no_info()
        self.controller._on_guide_toggle({"new": True})
        self._assert_guide_is_on()
        self.controller._on_about_toggle({"new": True})
        self._assert_about_is_on()
        self.controller._on_guide_toggle({"new": True})
        self._assert_guide_is_on()
        self.controller._on_guide_toggle({"new": False})
        self._assert_no_info()

    def _assert_guide_is_on(self):
        """Assert guide is on."""
        assert len(self.view.info_container.children) == 2
        assert self.view.guide in self.view.info_container.children

    def _assert_about_is_on(self):
        """Assert about is on."""
        assert len(self.view.info_container.children) == 1
        assert self.view.about in self.view.info_container.children

    def _assert_no_info(self):
        """Assert no info is shown."""
        assert len(self.view.info_container.children) == 0

    def _instansiate_mvc_components(self):
        """Instansiate `AppWrapper` MVC components."""
        self.model = AppWrapperModel()
        self.view = AppWrapperView()
        self.controller = AppWrapperContoller(self.model, self.view)
