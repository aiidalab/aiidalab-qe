from aiidalab_qe.app import AppController, AppModel, AppView
from aiidalab_qe.app.structure import StructureStepModel
from aiidalab_qe.common.wizard import State


class TestApp:
    @classmethod
    def setup_class(cls):
        cls.model = AppModel()
        cls.view = AppView()
        cls.controller = AppController(cls.model, cls.view)

    def test_load_app_from_process(self, generate_qeapp_workchain):
        """Test loading app from process."""
        workchain = generate_qeapp_workchain(
            relax_type="positions",
            spin_type="collinear",
            run_bands=True,
            run_pdos=False,
            functional="PBE",
        )
        self.model.process_uuid = workchain.node.uuid
        self.controller.load_wizard()
        wizard = self.controller.wizard
        assert wizard.configure_model.relax_type == "positions"
        assert wizard.configure_model.get_model("workchain").spin_type == "collinear"
        assert wizard.configure_model.get_model("bands").include is True
        assert wizard.configure_model.get_model("pdos").include is False
        assert wizard.configure_model.state == State.SUCCESS
        advanced_model = wizard.configure_model.get_model("advanced")
        pseudos_model = advanced_model.get_model("pseudos")
        assert len(pseudos_model.dictionary) > 0
        assert pseudos_model.functional == "PBE"

    def test_enable_toggles(self):
        """Test enable_toggles method."""
        assert self.view.guide_toggle.disabled is True
        assert self.view.about_toggle.disabled is True
        self.controller.enable_toggles()
        assert self.view.guide_toggle.disabled is False
        assert self.view.about_toggle.disabled is False

    def test_guide_toggle(self):
        """Test guide_toggle method."""
        self.controller.enable_toggles()
        self.controller._on_guide_toggle({"new": True})
        self._assert_guide_is_on()
        self.controller._on_guide_toggle({"new": False})
        self._assert_no_info()

    def test_about_toggle(self):
        """Test about_toggle method."""
        self.controller.enable_toggles()
        self.controller._on_about_toggle({"new": True})
        self._assert_about_is_on()
        self.controller._on_about_toggle({"new": False})
        self._assert_no_info()

    def test_toggle_switch(self):
        """Test toggle_switch method."""
        self.controller.enable_toggles()
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


def test_selecting_new_structure_unconfirms_model(generate_structure_data):
    model = StructureStepModel()
    model.structure_uuid = generate_structure_data().uuid
    assert model.has_structure
    model.confirm()
    model.structure_uuid = generate_structure_data().uuid
    assert not model.confirmed
