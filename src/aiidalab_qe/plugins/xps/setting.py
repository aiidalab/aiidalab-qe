# -*- coding: utf-8 -*-
"""Panel for Bands plugin.

"""
import ipywidgets as ipw
import traitlets as tl
from aiida.orm import Group, QueryBuilder, StructureData

from aiidalab_qe.common.panel import Panel


def install_pseudos(pseudo_group="pseudo_demo_pbe"):
    import os
    from pathlib import Path
    from subprocess import run

    base_url = "https://github.com/superstar54/xps-data/raw/main/pseudo_demo/"
    url = base_url + pseudo_group + ".aiida"

    env = os.environ.copy()
    env["PATH"] = f"{env['PATH']}:{Path.home().joinpath('.local', 'bin')}"

    def run_(*args, **kwargs):
        return run(*args, env=env, capture_output=True, check=True, **kwargs)

    run_(["verdi", "archive", "import", url])


class Setting(Panel):
    title = "XPS Settings"
    identifier = "xps"
    input_structure = tl.Instance(StructureData, allow_none=True)
    protocol = tl.Unicode(allow_none=True)

    core_hole_treatment_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Core hole treatment</h4></div>"""
    )
    core_hole_treatment_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
        You have three options:<br>
        (1) XCH(smear): places the excited electron into the conduction band, suitable for extend system.<br>
        (2) XCH(fixed): places the excited electron into the conduction band, suitable for extend system.<br>
        (3) Full: remove one electron from the system, suitable for molecule. </div>"""
    )

    pseudo_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Pseudopotential</h4></div>"""
    )
    pseudo_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
        Please select a pseudopotential group. Ground-state and excited-state pseudopotentials for each absorbing element.
        </div>"""
    )

    element_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Select element</h4></div>"""
    )
    element_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
        The list of elements to be considered for analysis. If no elements list is given, we instead calculate all elements in the structure.
        </div>"""
    )
    structure_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Structure</h4></div>"""
    )
    structure_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Below you can indicate if the material should be treated as a molecule
        or a crystal.
        </div>"""
    )
    supercell_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Cell size</h4></div>"""
    )
    supercell_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Define the minimum cell length in angstrom for the resulting supercell, and thus all output
        structures. The default value of 8.0 angstrom will be used
        if no input is given. Setting this value to 0.0 will
        instruct the CF to not scale up the input structure.
        </div>"""
    )
    binding_energy_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Absolute binding energy</h4></div>"""
    )
    binding_energy_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        To calculate the absolute binding energy, you need to provide the correction energy for the core electrons. The correction energy is Ecorr = E_core_hole - E_gipaw, where E_core_hole and E_gipaw are calculated by Etot - Etotps. Etot and Etotps can be found in the output when generating the pseudo potential. A offset corretion by fitting the experimental data is also added. Here is a example: C:339.79,O:668.22,F:955.73,Si:153.19
        </div>"""
    )

    def __init__(self, **kwargs):
        # Core hole treatment type
        self.core_hole_treatment = ipw.ToggleButtons(
            options=[
                ("XCH(smear)", "xch_smear"),
                ("XCH(fixed)", "xch_fixed"),
                ("Full", "full"),
            ],
            value="xch_smear",
        )
        self.pseudo_group = ipw.Dropdown(
            options=["pseudo_demo_pbe", "pseudo_demo_pbesol"],
            value="pseudo_demo_pbe",
            description="Group:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.elements_list = ipw.VBox()

        self.structure_type = ipw.ToggleButtons(
            options=[
                ("Molecule", "molecule"),
                ("Crystal", "crystal"),
            ],
            value="crystal",
        )
        self.supercell_min_parameter = ipw.FloatText(
            value=8.0,
            description="The minimum cell length (Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.calc_binding_energy = ipw.Checkbox(
            description="Calculate binding energy: ",
            indent=False,
            value=False,
        )

        self.children = [
            self.structure_title,
            self.structure_help,
            ipw.HBox(
                [self.structure_type],
            ),
            # self.core_hole_treatment_title,
            # self.core_hole_treatment_help,
            # self.core_hole_treatment,
            self.pseudo_title,
            self.pseudo_help,
            self.pseudo_group,
            self.element_title,
            self.element_help,
            ipw.HBox(
                [self.elements_list],
            ),
            # self.supercell_title,
            # self.supercell_help,
            # ipw.HBox(
            # [self.supercell_min_parameter],
            # ),
        ]
        self.pseudo_group.observe(self._update_pseudo, names="value")
        super().__init__(**kwargs)

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        elements_list = [
            element.description
            for element in self.elements_list.children
            if element.value
        ]
        if len(elements_list) == 0:
            elements_list = [
                element.description for element in self.elements_list.children
            ]
        if len(elements_list) == 0:
            raise Exception(
                f"No element is supported by pseudo group {self.pseudo_group.value}."
            )
        parameters = {
            "core_hole_treatment": self.core_hole_treatment.value,
            "structure_type": self.structure_type.value,
            "pseudo_group": self.pseudo_group.value,
            "correction_energies": self.correction_energies,
            "elements_list": elements_list,
        }
        return parameters

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.pseudo_group.value = input_dict.get("pseudo_group", "pseudo_demo_pbe")
        self.core_hole_treatment.value = input_dict.get(
            "core_hole_treatment", "xch_smear"
        )
        self.structure_type.value = input_dict.get("structure_type", "crystal")
        elements_list = input_dict.get("elements_list", {})
        for ele in self.elements_list.children:
            if ele.value in elements_list:
                ele.value = True

    @tl.observe("input_structure")
    def _update_structure(self, _=None):
        self._update_element_list()

    def _update_element_list(self):
        if self.input_structure is None:
            return
        structure = self.input_structure
        elements_list = [Kind.symbol for Kind in structure.kinds]
        checkbox_list = []
        qb = QueryBuilder()
        qb.append(Group, filters={"label": self.pseudo_group.value})
        if len(qb.all()) == 0:
            install_pseudos(self.pseudo_group.value)
        group = qb.all()[0][0]
        self.correction_energies = group.base.extras.get("correction")
        supported_elements = {}
        for key in self.correction_energies:
            ele, orbital = key.split("_")
            if ele not in supported_elements:
                supported_elements[ele] = [key]
            else:
                supported_elements[ele].append(key)
        # print("supported_elements: ", supported_elements)
        for ele in elements_list:
            if ele in supported_elements:
                for orbital in supported_elements[ele]:
                    checkbox_list += (
                        ipw.Checkbox(
                            description=orbital,
                            indent=False,
                            value=False,
                            layout=ipw.Layout(max_width="100%"),
                        ),
                    )
            else:
                checkbox_list += (
                    ipw.Checkbox(
                        description=f"{ele}, not supported by the selected pseudo group",
                        indent=False,
                        value=False,
                        disabled=True,
                        style={"description_width": "initial"},
                        layout=ipw.Layout(max_width="100%"),
                    ),
                )
        self.elements_list.children = checkbox_list

    def _update_pseudo(self, change):
        pseudo_group = change["new"]
        qb = QueryBuilder()
        qb.append(Group, filters={"label": pseudo_group})
        if len(qb.all()) == 0:
            install_pseudos(pseudo_group)
        self._update_element_list()
