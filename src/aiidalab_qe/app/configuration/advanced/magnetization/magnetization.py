import ipywidgets as ipw

from ..subsettings import AdvancedSubSettings
from .model import MagnetizationModel


class MagnetizationSettings(AdvancedSubSettings):
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

    identifier = "magnetization"

    def __init__(self, model: MagnetizationModel, **kwargs):
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

        self.description = ipw.HTML("<b>Magnetization:</b>")

        self.magnetization_type_toggle = ipw.ToggleButtons(
            options=[
                ("Starting Magnetization", "starting_magnetization"),
                ("Tot. Magnetization", "tot_magnetization"),
            ],
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "type"),
            (self.magnetization_type_toggle, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.magnetization_type_toggle, "disabled"),
            lambda override: not override,
        )

        self.tot_magnetization = ipw.BoundedIntText(
            min=0,
            max=100,
            step=1,
            disabled=True,
            description="Total magnetization:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "total"),
            (self.tot_magnetization, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.tot_magnetization, "disabled"),
            lambda override: not override,
        )

        self.elements_widget = ipw.VBox()

        self.container = ipw.VBox(
            children=[
                self.tot_magnetization,
            ]
        )

        self.children = []

        self.rendered = True

        self.refresh(which="all")

    def _on_input_structure_change(self, _):
        self.refresh(which="structure")

    def _on_electronic_type_change(self, _):
        self._switch_widgets()

    def _on_spin_type_change(self, _):
        self.refresh(which="spin")

    def _on_magnetization_type_change(self, _):
        self._toggle_widgets()

    def _update(self, which):
        if self.updated:
            return
        self._show_loading()
        self._model.update(which)
        self._build_kinds_widget()
        self._switch_widgets()
        self._toggle_widgets()
        self.updated = True

    def _show_loading(self):
        if self.rendered:
            self.elements_widget.children = [self.loading_message]

    def _build_kinds_widget(self):
        if not self.rendered:
            return

        children = []

        symbols = (
            self._model.input_structure.get_kind_names()
            if self._model.input_structure
            else []
        )

        for symbol in symbols:
            element_widget = ipw.BoundedFloatText(
                description=symbol,
                min=-4,
                max=4,
                step=0.1,
                disabled=True,
            )
            link = ipw.link(
                (self._model, "moments"),
                (element_widget, "value"),
                [
                    lambda d, symbol=symbol: d.get(symbol, 0.0),
                    lambda v, symbol=symbol: {
                        **self._model.moments,
                        symbol: v,
                    },
                ],
            )
            self.links.append(link)
            ipw.dlink(
                (self._model, "override"),
                (element_widget, "disabled"),
                lambda override: not override,
            )
            children.append(element_widget)

        self.elements_widget.children = children

    def _switch_widgets(self):
        if not self.rendered:
            return
        if self._model.spin_type == "none":
            children = []
        else:
            children = [self.description]
            if self._model.electronic_type == "metal":
                children.extend([self.magnetization_type_toggle, self.container])
            else:
                children.append(self.tot_magnetization)
        self.children = children

    def _toggle_widgets(self):
        if self._model.spin_type == "none" or not self.rendered:
            return
        self.container.children = [
            self.tot_magnetization
            if self._model.type == "tot_magnetization"
            else self.elements_widget
        ]
