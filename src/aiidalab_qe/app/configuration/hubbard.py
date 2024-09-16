import ipywidgets as ipw
import traitlets as tl
from pymatgen.core.periodic_table import Element

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.common.widgets import LoadingWidget

from .model import config_model as model


class HubbardSettings(ipw.VBox):
    """Widget for setting up Hubbard parameters."""

    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    def __init__(self, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading Hubbard settings widget")],
            **kwargs,
        )

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
        self._define_elements()
        self._build_hubbard_widget()
        if self._needs_eigenvalues_widget:
            self._define_default_eigenvalues()
            self._build_eigenvalues_widget()
        else:
            self._unsubscribe_eigenvalues_widget()
            self.eigenvalues_widget.children = []
        if isinstance(change["new"], HubbardStructureData):
            self._set_parameters_from_hubbard_structure()

    def _define_elements(self):
        if model.input_structure is None:
            self.elements = []
        else:
            self.elements = [
                *filter(
                    lambda element: (
                        element.is_transition_metal
                        or element.is_lanthanoid
                        or element.is_actinoid
                    ),
                    [Element(kind.symbol) for kind in model.input_structure.kinds],
                )
            ]
        self._needs_eigenvalues_widget = len(self.elements) > 0

    def _build_hubbard_widget(self):
        """Build the widget for defining Hubbard U values
        for each atomic species in the input structure.

        Returns:
            hubbard_widget (ipywidgets.VBox):
                The widget containing the input fields for defining Hubbard U values.
        """

        def get_manifold(element):
            """Get the Hubbard manifold for a given element.

            Parameters:
                element (Element):
                    The element for which to determine the Hubbard manifold.

            Returns:
                hubbard_manifold (str):
                    The Hubbard manifold for the given element.
            """
            valence = [
                orbital
                for orbital in element.electronic_structure.split(".")
                if "[" not in orbital
            ]
            orbital_shells = [shell[:2] for shell in valence]

            def is_condition_met(shell):
                return condition and condition in shell

            # Conditions for determining the Hubbard manifold
            # to be selected from the electronic structure
            conditions = {
                element.is_transition_metal: "d",
                element.is_lanthanoid or element.is_actinoid: "f",
                element.is_post_transition_metal
                or element.is_metalloid
                or element.is_halogen
                or element.is_chalcogen
                or element.symbol in ["C", "N", "P"]: "p",
                element.is_alkaline or element.is_alkali or element.is_noble_gas: "s",
            }

            condition = next(
                (shell for condition, shell in conditions.items() if condition), None
            )

            hubbard_manifold = next(
                (shell for shell in orbital_shells if is_condition_met(shell)), None
            )

            return hubbard_manifold

        def get_labels():
            """Get a list of labels for the Hubbard widget.

            Returns:
                labels (list):
                    A list of labels in the format "{kind} - {manifold}".
            """
            kind_list = model.input_structure.get_kind_names()
            hubbard_manifold_list = [
                get_manifold(Element(kind.symbol))
                for kind in model.input_structure.kinds
            ]
            labels = [
                f"{kind} - {manifold}"
                for kind, manifold in zip(kind_list, hubbard_manifold_list)
            ]
            return labels

        children = []
        self.hubbard_widget_links = []

        if model.input_structure is None:
            input_labels = []
        else:
            input_labels = get_labels()
            children.append(ipw.HTML("Define U value [eV] "))

        for label in input_labels:
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

        if self._needs_eigenvalues_widget:
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

    def _define_default_eigenvalues(self):
        model.hubbard.eigenvalues = [
            [
                [
                    [state + 1, spin, element.symbol, "-1"]  # default eigenvalue
                    for state in range(5 if element.is_transition_metal else 7)
                ]
                for spin in range(2)  # spin up and down
            ]
            for element in self.elements  # transition metals and lanthanoids
        ]

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

        children = []
        self.eigenvalues_widget_links = []

        for ei, element in enumerate(self.elements):
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

    def _set_parameters_from_hubbard_structure(self):
        hubbard_parameters = model.input_structure.hubbard.dict()["parameters"]
        parameters = {
            f"{model.input_structure.sites[item['atom_index']].kind_name} - {item['atom_manifold']}": item[
                "value"
            ]
            for item in hubbard_parameters
        }
        model.hubbard.parameters = parameters
        model.hubbard.activate = True

    def _on_hubbard_check(self, change):
        self.container.children = [self.hubbard_widget] if change["new"] else []

    def _on_eigenvalues_check(self, change):
        self.hubbard_widget.children = (
            [
                *self.hubbard_widget.children,
                self.eigenvalues_widget,
            ]
            if change["new"]
            else [*self.hubbard_widget.children][:-1]
        )
