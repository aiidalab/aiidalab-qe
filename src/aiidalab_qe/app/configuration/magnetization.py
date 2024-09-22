import ipywidgets as ipw

from .model import MagnetizationModel


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

    def __init__(self, model: MagnetizationModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[LoadingWidget("Loading magnetization settings widget")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_electronic_type_change,
            "electronic_type",
        )
        self._model.observe(
            self._on_magnetization_type_change,
            "magnetization_type",
        )

        self.kind_widget_links = []

        self.rendered = False

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

        self.kinds = ipw.VBox()

        self.container = ipw.VBox(
            children=[
                self.tot_magnetization,
            ]
        )

        self.children = [self.description]

        self.rendered = True

    def reset(self):
        self._model.reset()

    def _on_input_structure_change(self, change):
        self._model.set_defaults_from_structure()
        self._build_kinds_widget(change)

    def _on_electronic_type_change(self, change):
        self._switch_widgets(change)

    def _on_magnetization_type_change(self, change):
        self._toggle_widgets(change)

    def _build_kinds_widget(self, change):
        children = []

        if (input_structure := change["new"]) is None:
            labels = []
            for link in self.kind_widget_links:
                link.unlink()
            self.kind_widget_links.clear()
        else:
            labels = input_structure.get_kind_names()

        for label in labels:
            kind_widget = ipw.BoundedFloatText(
                description=label,
                min=-4,
                max=4,
                step=0.1,
                disabled=True,
            )
            link = ipw.link(
                (self._model, "moments"),
                (kind_widget, "value"),
                [
                    lambda d, label=label: d.get(label, 0.0),
                    lambda v, label=label: {
                        **self._model.moments,
                        label: v,
                    },
                ],
            )
            self.kind_widget_links.append(link)
            ipw.dlink(
                (self._model, "override"),
                (kind_widget, "disabled"),
                lambda override: not override,
            )
            children.append(kind_widget)

        self.kinds.children = children

    def _switch_widgets(self, change):
        children = [self.description]
        if change["new"] == "metal":
            children.extend([self.magnetization_type_toggle, self.container])
        else:
            children.append(self.tot_magnetization)
        self.children = children

    def _toggle_widgets(self, change):
        if change["new"] == "tot_magnetization":
            self.container.children = [self.tot_magnetization]
        else:
            self.container.children = [self.kinds]
