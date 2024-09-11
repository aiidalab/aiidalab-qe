import os

import ipywidgets as ipw
import traitlets as tl
from IPython.display import clear_output, display

from aiida import orm


class MagnetizationSettings(ipw.VBox):
    """Widget to set the type of magnetization used in the calculation:
    1) Tot_magnetization: Total majority spin charge - minority spin charge.
    2) Starting magnetization: Starting spin polarization on atomic type 'i' in a spin polarized (LSDA or noncollinear/spin-orbit) calculation.

    For Starting magnetization you can set each kind names defined in the StructureData (StructureDtaa.get_kind_names())
    Usually these are the names of the elements in the StructureData
    (For example 'C' , 'N' , 'Fe' . However the StructureData can have defined kinds like 'Fe1' and 'Fe2')
    The widget generate a dictionary that can be used to set initial_magnetic_moments in the builder of PwBaseWorkChain

    Attributes:
        input_structure(StructureData): trait that containes the input_strucgure (confirmed structure from previous step)
    """

    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    electronic_type = tl.Unicode()
    disabled = tl.Bool()
    _DEFAULT_TOT_MAGNETIZATION = 0.0
    _DEFAULT_DESCRIPTION = "<b>Magnetization: Input structure not confirmed</b>"

    def __init__(self, **kwargs):
        self.input_structure = orm.StructureData()
        self.input_structure_labels = []
        self.tot_magnetization = ipw.BoundedIntText(
            min=0,
            max=100,
            step=1,
            value=self._DEFAULT_TOT_MAGNETIZATION,
            disabled=True,
            description="Total magnetization:",
            style={"description_width": "initial"},
        )
        self.magnetization_type = ipw.ToggleButtons(
            options=[
                ("Starting Magnetization", "starting_magnetization"),
                ("Tot. Magnetization", "tot_magnetization"),
            ],
            value="starting_magnetization",
            style={"description_width": "initial"},
        )
        self.description = ipw.HTML(self._DEFAULT_DESCRIPTION)
        self.kinds = self.create_kinds_widget()
        self.kinds_widget_out = ipw.Output()
        self.magnetization_out = ipw.Output()
        self.magnetization_type.observe(self._render, "value")
        super().__init__(
            children=[
                self.description,
                self.magnetization_out,
                self.kinds_widget_out,
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    @tl.observe("disabled")
    def _disabled_changed(self, _):
        """Disable the widget"""
        if hasattr(self.kinds, "children") and self.kinds.children:
            for i in range(len(self.kinds.children)):
                self.kinds.children[i].disabled = self.disabled
        self.tot_magnetization.disabled = self.disabled
        self.magnetization_type.disabled = self.disabled

    def reset(self):
        self.disabled = True
        self.tot_magnetization.value = self._DEFAULT_TOT_MAGNETIZATION
        #
        if self.input_structure is None:
            self.description.value = self._DEFAULT_DESCRIPTION
            self.kinds = None
        else:
            self.description.value = "<b>Magnetization</b>"
            self.kinds = self.create_kinds_widget()

    def create_kinds_widget(self):
        if self.input_structure_labels:
            widgets_list = []
            for kind_label in self.input_structure_labels:
                kind_widget = ipw.BoundedFloatText(
                    description=kind_label,
                    min=-4,
                    max=4,
                    step=0.1,
                    value=0.0,
                    disabled=True,
                )
                widgets_list.append(kind_widget)
            kinds_widget = ipw.VBox(widgets_list)
        else:
            kinds_widget = None

        return kinds_widget

    @tl.observe("electronic_type")
    def _electronic_type_changed(self, change):
        with self.magnetization_out:
            clear_output()
            if change["new"] == "metal":
                display(self.magnetization_type)
                self._render({"new": self.magnetization_type.value})
            else:
                display(self.tot_magnetization)
                with self.kinds_widget_out:
                    clear_output()

    def update_kinds_widget(self):
        self.input_structure_labels = self.input_structure.get_kind_names()
        self.kinds = self.create_kinds_widget()
        self.description.value = "<b>Magnetization</b>"

    def _render(self, value):
        if value["new"] == "tot_magnetization":
            with self.kinds_widget_out:
                clear_output()
                display(self.tot_magnetization)
        else:
            self.display_kinds()

    def display_kinds(self):
        if "PYTEST_CURRENT_TEST" not in os.environ and self.kinds:
            with self.kinds_widget_out:
                clear_output()
                display(self.kinds)

    def _update_widget(self, change):
        self.input_structure = change["new"]
        self.update_kinds_widget()
        self.display_kinds()

    def get_magnetization(self):
        """Method to generate the dictionary with the initial magnetic moments"""
        magnetization = {}
        for i in range(len(self.kinds.children)):
            magnetization[self.input_structure_labels[i]] = self.kinds.children[i].value
        return magnetization

    def _set_magnetization_values(self, magnetic_moments):
        """Set magnetization"""
        # self.override.value = True
        with self.hold_trait_notifications():
            for i in range(len(self.kinds.children)):
                if isinstance(magnetic_moments, dict):
                    self.kinds.children[i].value = magnetic_moments.get(
                        self.kinds.children[i].description, 0.0
                    )
                else:
                    self.kinds.children[i].value = magnetic_moments

    def _set_tot_magnetization(self, tot_magnetization):
        """Set the total magnetization"""
        self.tot_magnetization.value = tot_magnetization
