import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

from .model import config_model as model


class HubbardSettings(ipw.VBox):
    """Widget for setting up Hubbard parameters."""

    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    def __init__(self, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading Hubbard settings widget")],
            **kwargs,
        )

        self.hubbard_widget_links = []
        self.eigenvalues_widget_links = []

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.activate_hubbard = ipw.Checkbox(
            description="",
            tooltip="Use Hubbard DFT+U.",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (model.hubbard, "activate"),
            (self.activate_hubbard, "value"),
        )
        ipw.dlink(
            (model, "override"),
            (self.activate_hubbard, "disabled"),
            lambda override: not override,
        )
        self.activate_hubbard.observe(self._on_hubbard_check, "value")

        self.eigenvalues_help = ipw.HTML(
            value="For transition metals and lanthanoids, the starting eigenvalues can be defined (Magnetic calculation).",
            layout=ipw.Layout(width="auto"),
        )
        self.eigenvalues_label = ipw.Checkbox(
            description="Define eigenvalues",
            tooltip="Define eigenvalues",
            indent=False,
            layout=ipw.Layout(max_width="30%"),
        )
        ipw.link(
            (model.hubbard, "eigenvalues_label"),
            (self.eigenvalues_label, "value"),
        )
        self.eigenvalues_label.observe(self._on_eigenvalues_check, "value")

        self.hubbard_widget = ipw.VBox()
        self.eigenvalues_widget = ipw.VBox()
        self.container = ipw.VBox()

        self.children = [
            ipw.HBox(
                children=[
                    ipw.HTML("<b>Hubbard (DFT+U)</b>"),
                    self.activate_hubbard,
                ]
            ),
            self.container,
        ]

        ipw.dlink(
            (model, "input_structure"),
            (self, "input_structure"),
        )

        self.rendered = True

    def reset(self):
        """Reset the widget."""
        self._unsubscribe_eigenvalues_widget()
        model.hubbard.reset()

    @tl.observe("input_structure")
    def _on_input_structure_change(self, change):
        self._unsubscribe_hubbard_widget()
        model.hubbard.set_defaults_from_structure(change["new"])
        self._build_hubbard_widget()
        if model.hubbard.needs_eigenvalues_widget:
            self._build_eigenvalues_widget()
        else:
            self._unsubscribe_eigenvalues_widget()
            self.eigenvalues_widget.children = []
        if isinstance(change["new"], HubbardStructureData):
            model.hubbard.set_parameters_from_hubbard_structure(change["new"])

    def _on_hubbard_check(self, change):
        self._toggle_hubbard_widget(change)

    def _on_eigenvalues_check(self, change):
        self._toggle_eigenvalues_widget(change)

    def _build_hubbard_widget(self):
        """Build the widget for defining Hubbard U values
        for each atomic species in the input structure.

        Returns:
            hubbard_widget (ipywidgets.VBox):
                The widget containing the input fields for defining Hubbard U values.
        """

        children = []

        if model.input_structure:
            children.append(ipw.HTML("Define U value [eV] "))

        for label in model.hubbard.input_labels:
            float_widget = ipw.BoundedFloatText(
                description=label,
                min=0,
                max=20,
                step=0.1,
                layout={"width": "160px"},
            )
            link = ipw.link(
                (model.hubbard, "parameters"),
                (float_widget, "value"),
                [
                    lambda p, label=label: p.get(label, 0.0),
                    lambda v, label=label: {**model.hubbard.parameters, label: v},
                ],
            )
            self.hubbard_widget_links.append(link)
            children.append(float_widget)

        if model.hubbard.needs_eigenvalues_widget:
            children.extend(
                [
                    self.eigenvalues_help,
                    self.eigenvalues_label,
                ]
            )

        self.hubbard_widget.children = children

    def _unsubscribe_hubbard_widget(self):
        for link in self.hubbard_widget_links:
            link.unlink()
        self.hubbard_widget_links.clear()

    def _build_eigenvalues_widget(self):
        """Build the widget for selecting eigenvalues of different kinds of atoms.

        Returns:
            ipywidgets.VBox: Widget for selecting eigenvalues.
        """

        def update(index, spin, state, symbol, value):
            """Update the eigenvalues list."""
            eigenvalues = [*model.hubbard.eigenvalues]
            eigenvalues[index][spin][state] = [state + 1, spin, symbol, value]
            return eigenvalues

        self._unsubscribe_eigenvalues_widget()

        children = []

        for ei, element in enumerate(model.hubbard.elements):
            es = element.symbol
            num_states = 5 if element.is_transition_metal else 7  # d or f states

            label_layout = ipw.Layout(justify_content="flex-start", width="50px")
            spin_up_row = ipw.HBox([ipw.Label("Up:", layout=label_layout)])
            spin_down_row = ipw.HBox([ipw.Label("Down:", layout=label_layout)])

            for si in range(num_states):
                eigenvalues_up = ipw.Dropdown(
                    description=f"{si+1}",
                    options=["-1", "0", "1"],
                    layout=ipw.Layout(width="65px"),
                    style={"description_width": "initial"},
                )
                link = ipw.link(
                    (model.hubbard, "eigenvalues"),
                    (eigenvalues_up, "value"),
                    [
                        lambda evs, ei=ei, si=si: evs[ei][0][si][-1],
                        lambda v, ei=ei, si=si, es=es: update(ei, 0, si, es, v),
                    ],
                )
                self.eigenvalues_widget_links.append(link)
                spin_up_row.children += (eigenvalues_up,)

                eigenvalues_down = ipw.Dropdown(
                    description=f"{si+1}",
                    options=["-1", "0", "1"],
                    layout=ipw.Layout(width="65px"),
                    style={"description_width": "initial"},
                )
                link = ipw.link(
                    (model.hubbard, "eigenvalues"),
                    (eigenvalues_down, "value"),
                    [
                        lambda evs, ei=ei, si=si: evs[ei][1][si][-1],
                        lambda v, ei=ei, si=si, es=es: update(ei, 1, si, es, v),
                    ],
                )
                self.eigenvalues_widget_links.append(link)
                spin_down_row.children += (eigenvalues_down,)

            children.append(
                ipw.HBox(
                    [
                        ipw.Label(element.symbol, layout=label_layout),
                        ipw.VBox(
                            children=[
                                spin_up_row,
                                spin_down_row,
                            ]
                        ),
                    ]
                )
            )

        self.eigenvalues_widget.children = children

    def _unsubscribe_eigenvalues_widget(self):
        for link in self.eigenvalues_widget_links:
            link.unlink()
        self.eigenvalues_widget_links.clear()

    def _toggle_hubbard_widget(self, change):
        self.container.children = [self.hubbard_widget] if change["new"] else []

    def _toggle_eigenvalues_widget(self, change):
        self.hubbard_widget.children = (
            [
                *self.hubbard_widget.children,
                self.eigenvalues_widget,
            ]
            if change["new"]
            else [*self.hubbard_widget.children][:-1]
        )
