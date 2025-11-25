import ipywidgets as ipw

from aiidalab_qe.common.panel import ConfigurationSettingsPanel
from aiidalab_qe.common.widgets import HBoxWithUnits

from .model import MagnetizationConfigurationSettingsModel


class MagnetizationConfigurationSettingsPanel(
    ConfigurationSettingsPanel[MagnetizationConfigurationSettingsModel],
):
    """Widget to set the type of magnetization used in the calculation:
    1) Total magnetization: Total majority spin charge - minority spin charge.
    2) Magnetic moments: Starting spin polarization on atomic type 'i' in a spin polarized (LSDA or noncollinear/spin-orbit) calculation.

    For Magnetic moments you can set each kind names defined in the StructureData (StructureData.get_kind_names())
    Usually these are the names of the elements in the StructureData
    (For example 'C' , 'N' , 'Fe' . However the StructureData can have defined kinds like 'Fe1' and 'Fe2')
    The widget generate a dictionary that can be used to set initial_magnetic_moments in the builder of PwBaseWorkChain

    Attributes:
        input_structure(StructureData): trait that contains the input_structure (confirmed structure from previous step)
    """

    def __init__(self, model: MagnetizationConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_input_structure_change,
            "structure_uuid",
        )
        self._model.observe(
            self._on_electronic_type_change,
            "electronic_type",
        )
        self._model.observe(
            self._on_spin_type_change,
            "spin_type",
        )
        self._model.observe(
            self._on_pseudos_dictionary_change,
            "dictionary",
        )
        self._model.observe(
            self._on_magnetization_type_change,
            "type",
        )

    def render(self):
        if self.rendered:
            return

        self.unit = "µ<sub>B</sub>"

        self.insulator_help = ipw.HTML("""
            <div style="line-height: 1.4;">
                <b>Note:</b> defining the starting magnetic moments per atomic species
                is available only for metallic systems. To enable the feature, please
                set the <b>Electronic type</b> to <b>Metal</b> in the <b>Basic
                settings</b> tab.
            </div>
        """)

        self.magnetization_type_help = ipw.HTML()
        ipw.dlink(
            (self._model, "type_help"),
            (self.magnetization_type_help, "value"),
        )

        self.magnetization_type = ipw.ToggleButtons(
            style={"button_width": "initial"},
            layout=ipw.Layout(margin="0 0 5px"),
        )
        ipw.dlink(
            (self._model, "type_options"),
            (self.magnetization_type, "options"),
        )
        ipw.link(
            (self._model, "type"),
            (self.magnetization_type, "value"),
        )

        self.tot_magnetization = ipw.BoundedIntText(
            min=0,
            max=100,
            step=1,
            description="Total magnetization:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "total"),
            (self.tot_magnetization, "value"),
        )

        self.tot_magnetization_with_unit = HBoxWithUnits(
            self.tot_magnetization,
            self.unit,
        )

        self.moments_list = ipw.VBox()

        self.container = ipw.VBox(
            children=[
                self.tot_magnetization_with_unit,
            ]
        )

        self.children = []

        self.rendered = True

        self.refresh(specific="widgets")

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_electronic_type_change(self, _):
        self._switch_widgets()
        self._model.update_type_help()

    def _on_spin_type_change(self, _):
        self.refresh(specific="spin")

    def _on_pseudos_dictionary_change(self, _):
        self.refresh(specific="dictionary")

    def _on_magnetization_type_change(self, _):
        self._toggle_widgets()
        self._model.update_type_help()

    def _update(self):
        self._show_loading()
        self._build_moments_list()
        self._switch_widgets()
        self._toggle_widgets()

    def _show_loading(self):
        if self.rendered:
            self.moments_list.children = [self.loading_message]

    def _build_moments_list(self):
        if not self.rendered:
            return

        children = []

        kind_names = (
            self._model.input_structure.get_kind_names()
            if self._model.has_structure
            else []
        )

        for kind_name in kind_names:
            kind_moment = ipw.BoundedFloatText(
                description=kind_name,
                min=-7,
                max=7,
                step=0.1,
                style={"description_width": "150px"},
            )
            link = ipw.link(
                (self._model, "moments"),
                (kind_moment, "value"),
                [
                    lambda moments, kind_name=kind_name: moments.get(kind_name, 0.0),
                    lambda value, kind_name=kind_name: {
                        **self._model.moments,
                        kind_name: value,
                    },
                ],
            )
            self._links.append(link)
            children.append(HBoxWithUnits(kind_moment, "µ<sub>B</sub>"))

        self.moments_list.children = children

    def _switch_widgets(self):
        if not self.rendered:
            return
        if self._model.spin_type == "none":
            children = []
        else:
            if self._model.electronic_type == "metal":
                children = [
                    self.magnetization_type,
                    self.magnetization_type_help,
                    self.container,
                ]

            else:
                children = [
                    self.insulator_help,
                    self.magnetization_type_help,
                    self.tot_magnetization_with_unit,
                ]
        self.children = children

    def _toggle_widgets(self):
        if self._model.spin_type == "none" or not self.rendered:
            return
        self.container.children = [
            self.tot_magnetization_with_unit
            if self._model.type == "tot_magnetization"
            else self.moments_list
        ]
