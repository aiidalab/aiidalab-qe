import ipywidgets as ipw

from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.common.widgets import LoadingWidget

from .model import ConfigurationModel


class HubbardSettings(ipw.VBox):
    """Widget for setting up Hubbard parameters."""

    def __init__(self, model: ConfigurationModel, **kwargs):
        self.loading_message = LoadingWidget("Loading Hubbard settings")

        super().__init__(
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.advanced.hubbard.observe(
            self._on_hubbard_activation,
            "is_active",
        )
        self._model.advanced.hubbard.observe(
            self._on_eigenvalues_definition,
            "has_eigenvalues",
        )

        self.links = []
        self.eigenvalues_widget_links = []

        self.rendered = False
        self.updated = False

    def render(self):
        if self.rendered:
            return

        self.activate_hubbard_checkbox = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (self._model.advanced.hubbard, "is_active"),
            (self.activate_hubbard_checkbox, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
            (self.activate_hubbard_checkbox, "disabled"),
            lambda override: not override,
        )

        self.eigenvalues_help = ipw.HTML(
            value="For transition metals and lanthanoids, the starting eigenvalues can be defined (Magnetic calculation).",
            layout=ipw.Layout(width="auto"),
        )
        self.define_eigenvalues_checkbox = ipw.Checkbox(
            description="Define eigenvalues",
            indent=False,
            layout=ipw.Layout(max_width="30%"),
        )
        ipw.link(
            (self._model.advanced.hubbard, "has_eigenvalues"),
            (self.define_eigenvalues_checkbox, "value"),
        )

        self.hubbard_widget = ipw.VBox()
        self.eigenvalues_widget = ipw.VBox()

        self.container = ipw.VBox()

        self.children = [
            ipw.HBox(
                children=[
                    ipw.HTML("<b>Hubbard (DFT+U)</b>"),
                    self.activate_hubbard_checkbox,
                ]
            ),
            self.container,
        ]

        self.rendered = True

        self.update()

    def update(self):
        if not self.updated:
            self._update()
            self.updated = True

    def reset(self):
        if not self._model.input_structure:
            self._unsubscribe()
        self._model.advanced.hubbard.reset()
        self.updated = False

    def _on_input_structure_change(self, _):
        self._unsubscribe()
        self._model.advanced.hubbard.update()
        self._build_hubbard_widget(rebuild=True)
        if isinstance(self._model.input_structure, HubbardStructureData):
            self._model.advanced.hubbard.set_parameters_from_hubbard_structure()

    def _on_hubbard_activation(self, _):
        self._model.advanced.hubbard.update()
        self._build_hubbard_widget()
        self._toggle_hubbard_widget()

    def _on_eigenvalues_definition(self, _):
        self._toggle_eigenvalues_widget()

    def _update(self, rebuild=False):
        self._unsubscribe()
        self._show_loading()
        self._model.advanced.hubbard.update()
        self._build_hubbard_widget(rebuild=rebuild)
        self._toggle_hubbard_widget()
        self._toggle_eigenvalues_widget()

    def _show_loading(self):
        if self.rendered:
            self.hubbard_widget.children = [self.loading_message]

    def _build_hubbard_widget(self, rebuild=False):
        if not self.rendered or len(self.hubbard_widget.children) > 1 and not rebuild:
            return

        children = []

        if self._model.input_structure and self._model.advanced.hubbard.is_active:
            children.append(ipw.HTML("Define U value [eV] "))

        for label in self._model.advanced.hubbard.orbital_labels:
            float_widget = ipw.BoundedFloatText(
                description=label,
                min=0,
                max=20,
                step=0.1,
                layout={"width": "160px"},
            )
            link = ipw.link(
                (self._model.advanced.hubbard, "parameters"),
                (float_widget, "value"),
                [
                    lambda p, label=label: p.get(label, 0.0),
                    lambda v, label=label: {
                        **self._model.advanced.hubbard.parameters,
                        label: v,
                    },
                ],
            )
            self.links.append(link)
            children.append(float_widget)

        if self._model.advanced.hubbard.needs_eigenvalues_widget:
            children.extend(
                [
                    self.eigenvalues_help,
                    self.define_eigenvalues_checkbox,
                ]
            )

        self.hubbard_widget.children = children

        if self._model.advanced.hubbard.needs_eigenvalues_widget:
            self._build_eigenvalues_widget()
        else:
            self.eigenvalues_widget.children = []

    def _build_eigenvalues_widget(self):
        def update(index, spin, state, symbol, value):
            """Update the eigenvalues list."""
            eigenvalues = [*self._model.advanced.hubbard.eigenvalues]
            eigenvalues[index][spin][state] = [state + 1, spin, symbol, value]
            return eigenvalues

        children = []

        for ei, element in enumerate(self._model.advanced.hubbard.applicable_elements):
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
                    (self._model.advanced.hubbard, "eigenvalues"),
                    (eigenvalues_up, "value"),
                    [
                        lambda evs, ei=ei, si=si: str(evs[ei][0][si][-1]),
                        lambda v, ei=ei, si=si, es=es: update(ei, 0, si, es, float(v)),
                    ],
                )
                self.links.append(link)
                spin_up_row.children += (eigenvalues_up,)

                eigenvalues_down = ipw.Dropdown(
                    description=f"{si+1}",
                    options=["-1", "0", "1"],
                    layout=ipw.Layout(width="65px"),
                    style={"description_width": "initial"},
                )
                link = ipw.link(
                    (self._model.advanced.hubbard, "eigenvalues"),
                    (eigenvalues_down, "value"),
                    [
                        lambda evs, ei=ei, si=si: str(evs[ei][1][si][-1]),
                        lambda v, ei=ei, si=si, es=es: update(ei, 1, si, es, float(v)),
                    ],
                )
                self.links.append(link)
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

    def _toggle_hubbard_widget(self):
        if not self.rendered:
            return
        widget = [self.hubbard_widget] if self._model.advanced.hubbard.is_active else []
        self.container.children = widget

    def _toggle_eigenvalues_widget(self):
        if not self.rendered:
            return
        self.hubbard_widget.children = (
            [
                *self.hubbard_widget.children,
                self.eigenvalues_widget,
            ]
            if self._model.advanced.hubbard.has_eigenvalues
            else [*self.hubbard_widget.children][:-1]
        )

    def _unsubscribe(self):
        for link in self.links:
            link.unlink()
        self.links.clear()
