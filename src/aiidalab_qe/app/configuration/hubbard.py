import ipywidgets as ipw
from IPython.display import clear_output, display
from pymatgen.core.periodic_table import Element

from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


class HubbardWidget(ipw.VBox):
    """Widget for setting up Hubbard parameters."""

    def __init__(self, input_structure=None):
        self.input_structure = input_structure
        self.eigenvalues_help = ipw.HTML(
            value="For transition metals and lanthanoids, the starting eigenvalues can be defined (Magnetic calculation).",
            layout=ipw.Layout(width="auto"),
        )
        self.activate_hubbard = ipw.Checkbox(
            description="",
            tooltip="Use Hubbard DFT+U.",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="10%"),
        )
        self.eigenvalues_label = ipw.Checkbox(
            description="Define eigenvalues",
            tooltip="Define eigenvalues",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="30%"),
        )
        self.hubbard_widget = self.create_hubbard_widget()
        self.hubbard_widget_out = ipw.Output()
        self.eigen_values_widget = self.create_eigenvalues_widget()
        self.eigen_values_widget_out = ipw.Output()

        super().__init__(
            children=[
                ipw.HBox(
                    children=[
                        ipw.HTML("<b>Hubbard (DFT+U)</b>"),
                        self.activate_hubbard,
                    ]
                ),
                self.hubbard_widget_out,
                self.eigen_values_widget_out,
            ]
        )
        self.activate_hubbard.observe(self.toggle_hubbard_widgets, names="value")
        self.eigenvalues_label.observe(self.toggle_eigenvalues_widgets, names="value")

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
        if self.input_structure is None:
            self.input_labels = []
        else:
            self.input_labels = self._hubbard_widget_labels()
            condition = _display_checkbox(self.input_structure.get_symbols_set())
        widgets_list = []
        for label in self.input_labels:
            hbox_container = ipw.HBox()
            float_widget = ipw.BoundedFloatText(
                description=label,
                min=0,
                max=20,
                step=0.1,
                value=0.0,
                layout={"width": "160px"},
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
        kind_list = self.input_structure.get_kind_names()
        hubbard_manifold_list = [
            self._get_hubbard_manifold(Element(x.symbol))
            for x in self.input_structure.kinds
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
        for kind in self.input_kinds_eigenvalues:
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
                    eigenvalues_down = ipw.Dropdown(
                        description=f"{i+1}",
                        options=["-1", "0", "1"],
                        layout=ipw.Layout(width="65px"),
                        style={"description_width": "initial"},
                    )
                    widgets_list_up.append(eigenvalues_up)
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
        occup_kinds_widget = ipw.VBox(kind_list)
        return occup_kinds_widget

    def update_widgets(self, change):
        """
        Update the widgets based on the given change.
        """
        self.input_structure = change
        self.hubbard_widget = self.create_hubbard_widget()
        self.eigenvalues_label.value = False
        if self.activate_hubbard.value:
            with self.hubbard_widget_out:
                clear_output()
                display(self.hubbard_widget)

        self.eigen_values_widget = self.create_eigenvalues_widget()
        if self.eigenvalues_label.value:
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
        self.set_hubbard_widget(parameters)
        self.activate_hubbard.value = True

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

    def _get_hubbard_u(self) -> dict:
        """
        Get the Hubbard U values for each input label.

        Returns:
            dict: A dictionary containing the Hubbard U values for each input label.
                    The dictionary format is {'kind_name - hubbard_manifold': U_value}

        """
        hubbard_u = {}
        for index, label in enumerate(self.input_labels):
            value_hubbard = self.hubbard_widget.children[index + 1].children[0].value
            if value_hubbard != 0:
                hubbard_u[label] = (
                    self.hubbard_widget.children[index + 1].children[0].value
                )
        return hubbard_u

    def _get_starting_ns_eigenvalue(self) -> list:
        """
        Get the starting ns eigenvalues for transition metal and lanthanoid elements.

        Returns:
            list: A list of starting ns eigenvalues for each element.
                  Each element in the list is a list containing the following information:
                  - The eigenvalue
                  - The spin index
                  - The element symbol
                  - The value of the eigenvalue
        """
        starting_ns_eigenvalue = []
        for index, kind in enumerate(self.input_kinds_eigenvalues):
            if kind[2].is_transition_metal or kind[2].is_lanthanoid:
                if kind[2].is_transition_metal:
                    num_states = 5
                else:
                    num_states = 7
                for i in range(2):  # up and down
                    spin = (
                        self.eigen_values_widget.children[index].children[1].children[i]
                    )
                    for j in range(num_states):
                        value_eigenvalue = int(spin.children[j + 1].value)
                        if value_eigenvalue != -1:
                            starting_ns_eigenvalue.append(
                                [j + 1, i + 1, kind[1], value_eigenvalue]
                            )

        return starting_ns_eigenvalue

    def set_hubbard_widget(self, parameters):
        """
        Set the Hubbard widget based on the given parameters.

        Parameters:
            parameters (dict): A dictionary containing the Hubbard U values for each input label.
                               The dictionary format is {'kind_name - hubbard_manifold': U_value}
        """
        for index, label in enumerate(self.input_labels):
            if label in parameters:
                self.hubbard_widget.children[index + 1].children[0].value = parameters[
                    label
                ]

    def set_eigenvalues_widget(self, parameters):
        """
        Set the eigenvalues widget based on the given parameters.

        Parameters:
            parameters (list): A list of starting ns eigenvalues for each element.
                               Each element in the list is a list containing the following information:
                               - The eigenvalue
                               - The spin index
                               - The element symbol
                               - The value of the eigenvalue
        """
        for param in parameters:
            eigenvalue, spin_index, element_symbol, value = param
            for index, kind in enumerate(self.input_kinds_eigenvalues):
                if kind[1] == element_symbol:
                    spin = (
                        self.eigen_values_widget.children[index]
                        .children[1]
                        .children[spin_index - 1]
                    )
                    spin.children[eigenvalue].value = str(value)

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

    @property
    def hubbard_dict(self) -> dict:
        if self.activate_hubbard.value:
            hubbard_dict = {
                "hubbard_u": self._get_hubbard_u(),
            }
        else:
            hubbard_dict = {}
        return hubbard_dict

    @property
    def eigenvalues_dict(self) -> dict:
        if self.eigenvalues_label.value:
            eigenvalues_dict = {
                "starting_ns_eigenvalue": self._get_starting_ns_eigenvalue()
            }
        else:
            eigenvalues_dict = {}
        return eigenvalues_dict
