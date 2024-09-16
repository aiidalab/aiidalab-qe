import ipywidgets as ipw
import traitlets as tl

from aiida import orm

from .model import config_model as model


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

    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    electronic_type = tl.Unicode()
    magnetization_type = tl.Unicode()

    def __init__(self, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[LoadingWidget("Loading magnetization settings widget")],
            **kwargs,
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
            (model.magnetization, "type"),
            (self.magnetization_type_toggle, "value"),
        )
        ipw.dlink(
            (model, "override"),
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
            (model.magnetization, "total"),
            (self.tot_magnetization, "value"),
        )
        ipw.dlink(
            (model, "override"),
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

        ipw.dlink(
            (model, "input_structure"),
            (self, "input_structure"),
        )
        ipw.dlink(
            (model, "electronic_type"),
            (self, "electronic_type"),
        )
        ipw.dlink(
            (model.magnetization, "type"),
            (self, "magnetization_type"),
        )

        self.rendered = True

    def reset(self):
        model.magnetization.reset()

    @tl.observe("input_structure")
    def _on_input_structure_change(self, change):
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
                (model.magnetization, "moments"),
                (kind_widget, "value"),
                [
                    lambda d, label=label: d.get(label, 0.0),
                    lambda v, label=label: {
                        **model.magnetization.moments,
                        label: v,
                    },
                ],
            )
            self.kind_widget_links.append(link)
            ipw.dlink(
                (model, "override"),
                (kind_widget, "disabled"),
                lambda override: not override,
            )
            children.append(kind_widget)

        self.kinds.children = children

    @tl.observe("electronic_type")
    def _on_electronic_type_change(self, change):
        children = [self.description]
        if change["new"] == "metal":
            children.extend([self.magnetization_type_toggle, self.container])
        else:
            children.append(self.tot_magnetization)
        self.children = children

    @tl.observe("magnetization_type")
    def _on_magnetization_type_change(self, change):
        if change["new"] == "tot_magnetization":
            self.container.children = [self.tot_magnetization]
        else:
            self.container.children = [self.kinds]
