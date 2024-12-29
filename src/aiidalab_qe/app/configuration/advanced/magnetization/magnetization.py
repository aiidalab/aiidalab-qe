import ipywidgets as ipw

from ..subsettings import AdvancedConfigurationSubSettingsPanel
from .model import MagnetizationConfigurationSettingsModel


class MagnetizationConfigurationSettingsPanel(
    AdvancedConfigurationSubSettingsPanel[MagnetizationConfigurationSettingsModel],
):
    """Widget to set the type of magnetization used in the calculation:
    1) Tot_magnetization: Total majority spin charge - minority spin charge.
    2) Starting magnetization: Starting spin polarization on atomic type 'i' in a spin polarized (LSDA or noncollinear/spin-orbit) calculation.

    For Starting magnetization you can set each kind names defined in the StructureData (StructureData.get_kind_names())
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
            "input_structure",
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
            self._on_magnetization_type_change,
            "type",
        )

    def render(self):
        if self.rendered:
            return

        self.description = ipw.HTML("""
            <div style="margin-bottom: 5px;">
                <b>Magnetization:</b>
                <br>
                The default starting magnetization is computed as the theoretical
                magnetic moment scaled by the valence charge defined in the selected
                pseudopotential family.
            </div>
        """)

        self.magnetization_type = ipw.ToggleButtons(
            style={
                "description_width": "initial",
                "button_width": "initial",
            },
            layout=ipw.Layout(margin="0 0 10px 0"),
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
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "total"),
            (self.tot_magnetization, "value"),
        )

        self.kind_moment_widgets = ipw.VBox()

        self.container = ipw.VBox(
            children=[
                self.tot_magnetization,
            ]
        )

        self.children = []

        self.rendered = True

        self.refresh(specific="widgets")

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_electronic_type_change(self, _):
        self._switch_widgets()

    def _on_spin_type_change(self, _):
        self.refresh(specific="spin")

    def _on_magnetization_type_change(self, _):
        self._toggle_widgets()

    def _update(self, specific=""):
        if self.updated:
            return
        self._show_loading()
        if not self._model.loaded_from_process or (specific and specific != "widgets"):
            self._model.update(specific)
        self._build_kinds_widget()
        self._switch_widgets()
        self._toggle_widgets()
        self.updated = True

    def _show_loading(self):
        if self.rendered:
            self.kind_moment_widgets.children = [self.loading_message]

    def _build_kinds_widget(self):
        if not self.rendered:
            return

        children = []

        kind_names = (
            self._model.input_structure.get_kind_names()
            if self._model.input_structure
            else []
        )

        for kind_name in kind_names:
            kind_moment_widget = ipw.BoundedFloatText(
                description=kind_name,
                min=-4,
                max=4,
                step=0.1,
            )
            link = ipw.link(
                (self._model, "moments"),
                (kind_moment_widget, "value"),
                [
                    lambda moments, kind_name=kind_name: moments.get(kind_name, 0.0),
                    lambda value, kind_name=kind_name: {
                        **self._model.moments,
                        kind_name: value,
                    },
                ],
            )
            self.links.append(link)
            children.append(kind_moment_widget)

        self.kind_moment_widgets.children = children

    def _switch_widgets(self):
        if not self.rendered:
            return
        if self._model.spin_type == "none":
            children = []
        else:
            children = [self.description]
            if self._model.electronic_type == "metal":
                children.extend(
                    [
                        self.magnetization_type,
                        self.container,
                    ]
                )
            else:
                children.append(self.tot_magnetization)
        self.children = children

    def _toggle_widgets(self):
        if self._model.spin_type == "none" or not self.rendered:
            return
        self.container.children = [
            self.tot_magnetization
            if self._model.type == "tot_magnetization"
            else self.kind_moment_widgets
        ]
