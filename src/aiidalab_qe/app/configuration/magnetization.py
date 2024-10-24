import ipywidgets as ipw

from aiidalab_qe.common.widgets import LoadingWidget

from .model import ConfigurationModel


class MagnetizationSettings(ipw.VBox):
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

    def __init__(self, model: ConfigurationModel, **kwargs):
        self.loading_message = LoadingWidget("Loading magnetization settings")

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.workchain.observe(
            self._on_electronic_type_change,
            "electronic_type",
        )
        self._model.workchain.observe(
            self._on_spin_type_change,
            "spin_type",
        )
        self._model.advanced.magnetization.observe(
            self._on_magnetization_type_change,
            "type",
        )

        self.links = []

        self.rendered = False
        self.updated = False

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
            (self._model.advanced.magnetization, "type"),
            (self.magnetization_type_toggle, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
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
            (self._model.advanced.magnetization, "total"),
            (self.tot_magnetization, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
            (self.tot_magnetization, "disabled"),
            lambda override: not override,
        )

        self.kinds_widget = ipw.VBox()

        self.container = ipw.VBox(
            children=[
                self.tot_magnetization,
            ]
        )

        self.children = [self.description]

        self.rendered = True

        self.update()

    def update(self):
        if not self.updated:
            self._update()
            self.updated = True

    def reset(self):
        if not self._model.input_structure:
            self._unsubscribe()
        self._model.advanced.magnetization.reset()
        self.updated = False

    def _on_input_structure_change(self, _):
        self._model.advanced.magnetization.update()
        self._update(rebuild=True)

    def _on_electronic_type_change(self, _):
        self._switch_widgets()

    def _on_spin_type_change(self, _):
        self._update()

    def _on_magnetization_type_change(self, _):
        self._toggle_widgets()

    def _update(self, rebuild=False):
        self._unsubscribe()
        self._show_loading()
        self._model.advanced.magnetization.update()
        self._build_kinds_widget(rebuild=rebuild)
        self._switch_widgets()
        self._toggle_widgets()

    def _show_loading(self):
        if self.rendered:
            self.kinds_widget.children = [self.loading_message]

    def _build_kinds_widget(self, rebuild=False):
        if not self.rendered or len(self.kinds_widget.children) > 0 and not rebuild:
            return

        children = []

        if (
            self._model.workchain.spin_type == "none"
            or self._model.input_structure is None
        ):
            labels = []
        else:
            labels = self._model.input_structure.get_kind_names()

        for label in labels:
            kind_widget = ipw.BoundedFloatText(
                description=label,
                min=-4,
                max=4,
                step=0.1,
                disabled=True,
            )
            link = ipw.link(
                (self._model.advanced.magnetization, "moments"),
                (kind_widget, "value"),
                [
                    lambda d, label=label: d.get(label, 0.0),
                    lambda v, label=label: {
                        **self._model.advanced.magnetization.moments,
                        label: v,
                    },
                ],
            )
            self.links.append(link)
            ipw.dlink(
                (self._model.advanced, "override"),
                (kind_widget, "disabled"),
                lambda override: not override,
            )
            children.append(kind_widget)

        self.kinds_widget.children = children

    def _switch_widgets(self):
        if not self.rendered:
            return
        if self._model.workchain.spin_type == "none":
            children = []
        else:
            children = [self.description]
            if self._model.workchain.electronic_type == "metal":
                children.extend([self.magnetization_type_toggle, self.container])
            else:
                children.append(self.tot_magnetization)
        self.children = children

    def _toggle_widgets(self):
        if self._model.workchain.spin_type == "none" or not self.rendered:
            return
        if self._model.advanced.magnetization.type == "tot_magnetization":
            self.container.children = [self.tot_magnetization]
        else:
            self.container.children = [self.kinds_widget]

    def _unsubscribe(self):
        for link in self.links:
            link.unlink()
        self.links.clear()
