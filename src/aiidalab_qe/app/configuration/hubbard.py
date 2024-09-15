import ipywidgets as ipw
import traitlets as tl
from IPython.display import clear_output, display
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

        self.eigenvalues_help = ipw.HTML(
            value="For transition metals and lanthanoids, the starting eigenvalues can be defined (Magnetic calculation).",
            layout=ipw.Layout(width="auto"),
        )
        self.activate_hubbard = ipw.Checkbox(
            description="",
            tooltip="Use Hubbard DFT+U.",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (model, "activate_hubbard"),
            (self.activate_hubbard, "value"),
        )
        ipw.dlink(
            (model, "override"),
            (self.activate_hubbard, "disabled"),
            lambda override: not override,
        )
        self.activate_hubbard.observe(self.toggle_hubbard_widgets, "value")

        self.eigenvalues_label = ipw.Checkbox(
            description="Define eigenvalues",
            tooltip="Define eigenvalues",
            indent=False,
            layout=ipw.Layout(max_width="30%"),
        )
        ipw.link(
            (model, "eigenvalues_label"),
            (self.eigenvalues_label, "value"),
        )
        self.eigenvalues_label.observe(self.toggle_eigenvalues_widgets, "value")

        self.hubbard_widget_out = ipw.Output()
        self.eigen_values_widget_out = ipw.Output()

        self.children = [
            ipw.HBox(
                children=[
                    ipw.HTML("<b>Hubbard (DFT+U)</b>"),
                    self.activate_hubbard,
                ]
            ),
            self.hubbard_widget_out,
            self.eigen_values_widget_out,
        ]

        ipw.dlink(
            (model, "input_structure"),
            (self, "input_structure"),
        )

        self.rendered = True

    def create_hubbard_widget(self):
        """
        Creates a widget for defining Hubbard U values for each atomic species in the input structure.

        Returns:
            hubbard_widget (ipywidgets.VBox): The widget containing the input fields for defining Hubbard U values.
        """

        def _display_checkbox(symbols):
            return any(
                Element(symbol).is_transition_metal
                or Element(symbol).is_lanthanoid
                or Element(symbol).is_actinoid
                for symbol in symbols
            )

        condition = False
        if model.input_structure is None:
            model.input_labels = []
        else:
            model.input_labels = self._hubbard_widget_labels()
            condition = _display_checkbox(model.input_structure.get_symbols_set())
        widgets_list = []
        for label in model.input_labels:
            hbox_container = ipw.HBox()
            float_widget = ipw.BoundedFloatText(
                description=label,
                min=0,
                max=20,
                step=0.1,
                layout={"width": "160px"},
            )
            ipw.dlink(
                (model, "hubbard_dict"),
                (float_widget, "value"),
                lambda d, label=label: d.get(label, 0.0),
            )
            float_widget.observe(
                lambda change, label=label: self._update_hubbard_dict(
                    {**model.hubbard_dict, label: change["new"]}
                ),
                "value",
            )
            hbox_container.children = [float_widget]
            widgets_list.append(hbox_container)

        if condition:
            hubbard_widget = ipw.VBox(
                [
                    ipw.HTML("Define U value [eV] "),
                    *widgets_list,
                    self.eigenvalues_help,
                    self.eigenvalues_label,
                ]
            )
        else:
            hubbard_widget = ipw.VBox([ipw.HTML("Define U value [eV] "), *widgets_list])
        return hubbard_widget

    def _hubbard_widget_labels(self):
        """
        Returns a list of labels for the Hubbard widget.

        The labels are generated based on the kind names and the corresponding Hubbard manifolds
        of the input structure.

        Returns:
            list: A list of labels in the format "{kind} - {manifold}".
        """
        kind_list = model.input_structure.get_kind_names()
        hubbard_manifold_list = [
            self._get_hubbard_manifold(Element(x.symbol))
            for x in model.input_structure.kinds
        ]
        result = [
            f"{kind} - {manifold}"
            for kind, manifold in zip(kind_list, hubbard_manifold_list)
        ]
        return result

    def _get_hubbard_manifold(self, element):
        """
        Get the Hubbard manifold for a given element.

        Parameters:
        element (Element): The element for which to determine the Hubbard manifold.

        Returns:
        str: The Hubbard manifold for the given element.
        """
        valence = [
            orbital
            for orbital in element.electronic_structure.split(".")
            if "[" not in orbital
        ]
        orbital_shells = [shell[:2] for shell in valence]

        def is_condition_met(shell):
            return condition and condition in shell

        # Conditions for determining the Hubbard manifold to be selected from the electronic structure
        hubbard_conditions = {
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
            (shell for condition, shell in hubbard_conditions.items() if condition),
            None,
        )

        hubbard_manifold = next(
            (shell for shell in orbital_shells if is_condition_met(shell)), None
        )

        return hubbard_manifold

    def _update_hubbard_dict(self, hubbard_dict):
        model.hubbard_dict = hubbard_dict

    def create_eigenvalues_widget(self):
        """
        Creates and returns a widget for selecting eigenvalues of different kinds of atoms.

        Returns:
        occup_kinds_widget (ipywidgets.VBox): Widget for selecting eigenvalues.
        """

        if self.input_structure is None:
            self.input_kinds_eigenvalues = []
        else:
            list_of_kinds = [
                [index + 1, value.name, Element(value.symbol)]
                for index, value in enumerate(self.input_structure.kinds)
            ]
            self.input_kinds_eigenvalues = [
                x
                for x in list_of_kinds
                if x[2].is_transition_metal or x[2].is_lanthanoid
            ]

        kind_list = []
        for index, kind in enumerate(self.input_kinds_eigenvalues):
            if kind[2].is_transition_metal:
                num_states = 5  # d states
            if kind[2].is_lanthanoid:
                num_states = 7  # f states

            if kind[2].is_transition_metal or kind[2].is_lanthanoid:
                widgets_list_up = []
                widgets_list_down = []
                for i in range(num_states):
                    eigenvalues_up = ipw.Dropdown(
                        description=f"{i+1}",
                        options=["-1", "0", "1"],
                        layout=ipw.Layout(width="65px"),
                        style={"description_width": "initial"},
                    )
                    ipw.dlink(
                        (model, "eigenvalues_list"),
                        (eigenvalues_up, "value"),
                        lambda evl, index=index, i=i: evl[index][0][i][-1],
                    )
                    eigenvalues_up.observe(
                        lambda change, index=index, i=i: self._update_eigenvalues_list(
                            index, 0, i, change["new"]
                        ),
                        "value",
                    )
                    widgets_list_up.append(eigenvalues_up)

                    eigenvalues_down = ipw.Dropdown(
                        description=f"{i+1}",
                        options=["-1", "0", "1"],
                        layout=ipw.Layout(width="65px"),
                        style={"description_width": "initial"},
                    )
                    ipw.dlink(
                        (model, "eigenvalues_list"),
                        (eigenvalues_down, "value"),
                        lambda evl, index=index, i=i: evl[index][1][i][-1],
                    )
                    eigenvalues_down.observe(
                        lambda change, index=index, i=i: self._update_eigenvalues_list(
                            index, 1, i, change["new"]
                        ),
                        "value",
                    )
                    widgets_list_down.append(eigenvalues_down)

                row_up = ipw.HBox(
                    children=[
                        ipw.Label(
                            "Up:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="50px"
                            ),
                        ),
                        *widgets_list_up,
                    ],
                )

                row_down = ipw.HBox(
                    children=[
                        ipw.Label(
                            "Down:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="50px"
                            ),
                        ),
                        *widgets_list_down,
                    ],
                )

                eigenvalues_container = ipw.VBox(children=[row_up, row_down])

                kind_container = ipw.HBox(
                    children=[
                        ipw.Label(
                            kind[1],
                            layout=ipw.Layout(
                                justify_content="flex-start", width="50px"
                            ),
                        ),
                        eigenvalues_container,
                    ]
                )

                kind_list.append(kind_container)

        return ipw.VBox(kind_list)

    def _update_eigenvalues_list(self, index, spin, eigenvalue, value):
        eigenvalues_list = {**model.eigenvalues_list}
        eigenvalues_list[index][spin][eigenvalue] = value
        model.eigenvalues_list = eigenvalues_list

    @tl.observe("input_structure")
    def update_widgets(self, _):
        """
        Update the widgets based on the given change.
        """
        self.hubbard_widget = self.create_hubbard_widget()
        model.eigenvalues_label = False  # TODO permanently disabled?
        if model.activate_hubbard:
            with self.hubbard_widget_out:
                clear_output()
                display(self.hubbard_widget)

        self.eigen_values_widget = self.create_eigenvalues_widget()
        if model.eigenvalues_label:
            with self.eigen_values_widget_out:
                clear_output()
                display(self.eigen_values_widget)

        if isinstance(self.input_structure, HubbardStructureData):
            self.set_parameters_from_hubbardstructure(self.input_structure)

    def set_parameters_from_hubbardstructure(self, hubbard_structure):
        hubbard_parameters = hubbard_structure.hubbard.dict()["parameters"]
        parameters = {
            f"{hubbard_structure.sites[item['atom_index']].kind_name} - {item['atom_manifold']}": item[
                "value"
            ]
            for item in hubbard_parameters
        }
        model.hubbard_dict = parameters
        model.activate_hubbard = True

    def toggle_hubbard_widgets(self, change):
        """
        Toggle the visibility of the Hubbard widgets based on the value of check box.
        """
        if change["new"]:
            with self.hubbard_widget_out:
                clear_output()
                display(self.hubbard_widget)
            if self.eigenvalues_label.value:
                with self.eigen_values_widget_out:
                    clear_output()
                    display(self.eigen_values_widget)
        else:
            with self.hubbard_widget_out:
                clear_output()
            with self.eigen_values_widget_out:
                clear_output()

    def toggle_eigenvalues_widgets(self, change):
        """
        Toggle the visibility of eigenvalues widgets based on the value of the eigenvalues check box.
        """
        if change["new"]:
            with self.eigen_values_widget_out:
                clear_output()
                display(self.eigen_values_widget)
        else:
            with self.eigen_values_widget_out:
                clear_output()

    def reset(self):
        """Reset the widget."""
        self.activate_hubbard.value = False
        self.eigenvalues_label.value = False
        self.hubbard_widget = self.create_hubbard_widget()
        self.eigen_values_widget = self.create_eigenvalues_widget()
        with self.hubbard_widget_out:
            clear_output()
        with self.eigen_values_widget_out:
            clear_output()
        if isinstance(self.input_structure, HubbardStructureData):
            self.set_parameters_from_hubbardstructure(self.input_structure)
