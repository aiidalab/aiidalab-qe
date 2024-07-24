"""Panel for XAS plugin."""

import os
import tarfile
from importlib import resources
from pathlib import Path

import ipywidgets as ipw
import requests
import traitlets as tl
import yaml
from aiida import orm

from aiidalab_qe.common.panel import Panel
from aiidalab_qe.plugins import xas as xas_folder

PSEUDO_TOC = yaml.safe_load(resources.read_text(xas_folder, "pseudo_toc.yaml"))
pseudo_data_dict = PSEUDO_TOC["pseudos"]
xch_elements = PSEUDO_TOC["xas_xch_elements"]

base_url = "https://github.com/PNOGillespie/Core_Level_Spectra_Pseudos/raw/main"
head_path = f"{Path.home()}/.local/lib"
dir_header = "cls_pseudos"
functionals = ["pbe"]
core_wfc_dir = "core_wfc_data"
gipaw_dir = "gipaw_pseudos"
ch_pseudo_dir = "ch_pseudos/star1s"


def _load_or_import_nodes_from_filenames(in_dict, path, core_wfc_data=False):
    for filename in in_dict.values():
        try:
            orm.load_node(filename)
        except BaseException:
            if not core_wfc_data:
                new_upf = orm.UpfData(f"{path}/{filename}", filename=filename)
                new_upf.label = filename
                new_upf.store()
            else:
                new_singlefile = orm.SinglefileData(
                    f"{path}/{filename}", filename="stdout"
                )
                new_singlefile.label = filename
                new_singlefile.store()


def _download_extract_pseudo_archive(func):
    target_dir = f"{head_path}/{dir_header}/{func}"
    archive_filename = f"{func}_ch_pseudos.tgz"
    remote_archive_filename = f"{base_url}/{func}/{archive_filename}"
    local_archive_filename = f"{target_dir}/{archive_filename}"

    env = os.environ.copy()
    env["PATH"] = f"{env['PATH']}:{Path.home() / '.local' / 'lib'}"

    response = requests.get(remote_archive_filename, timeout=30)
    response.raise_for_status()
    with open(local_archive_filename, "wb") as handle:
        handle.write(response.content)
        handle.flush()
        response.close()

    with tarfile.open(local_archive_filename, "r:gz") as tarfil:
        tarfil.extractall(target_dir)


url = f"{base_url}"
for func in functionals:
    target_dir = f"{head_path}/{dir_header}/{func}"
    os.makedirs(target_dir, exist_ok=True)
    archive_filename = f"{func}_ch_pseudos.tgz"
    archive_found = False
    for entry in os.listdir(target_dir):
        if entry == archive_filename:
            archive_found = True
    if not archive_found:
        _download_extract_pseudo_archive(func)


# Check all the pseudos/core-wfc data files in the TOC dictionary
# and load/check all of them before proceeding. Note that this
# approach relies on there not being multiple instances of nodes
# with the same label.
for func in functionals:
    gipaw_pseudo_dict = pseudo_data_dict[func]["gipaw_pseudos"]
    core_wfc_dict = pseudo_data_dict[func]["core_wavefunction_data"]
    core_hole_pseudo_dict = pseudo_data_dict[func]["core_hole_pseudos"]
    main_path = f"{head_path}/{dir_header}/{func}"
    core_wfc_dir = f"{main_path}/core_wfc_data"
    gipaw_dir = f"{main_path}/gipaw_pseudos"
    ch_pseudo_dir = f"{main_path}/ch_pseudos/star1s"
    # First, check that the local directories contain what's in the pseudo_toc
    for pseudo_dir, pseudo_dict in zip(
        [gipaw_dir, core_wfc_dir, ch_pseudo_dir],
        [gipaw_pseudo_dict, core_wfc_dict, core_hole_pseudo_dict],
    ):
        pseudo_toc_mismatch = os.listdir(pseudo_dir) != pseudo_dict.values()

    # Re-download the relevant archive if there is a mismatch
    if pseudo_toc_mismatch:
        _download_extract_pseudo_archive(func)

    _load_or_import_nodes_from_filenames(
        in_dict=gipaw_pseudo_dict,
        path=gipaw_dir,
    )
    _load_or_import_nodes_from_filenames(
        in_dict=core_wfc_dict, path=core_wfc_dir, core_wfc_data=True
    )
    _load_or_import_nodes_from_filenames(
        in_dict=core_hole_pseudo_dict["1s"], path=ch_pseudo_dir
    )


class Setting(Panel):
    title = "XAS Settings"
    identifier = "xas"
    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    protocol = tl.Unicode(allow_none=True)

    element_selection_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Element and Core-Hole Treatment Setting.</h4></div>"""
    )

    # TODO: The element selection should lock the "Confirm" button if no elements have been
    # selected for XAS calculation.

    element_selection_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
        To select elements for calculation of K-edge spectra:<br>
        (1) Tick the checkbox for each element symbol to select the element for calculation.<br>
        (2) Select the core-hole treatment scheme from the dropdown box.<br>
        <br>
        There are three supported options for core-hole treatment:<br>
         - FCH: Remove one electron from the system (any occupations scheme).<br>
         - XCH (Smearing): places the excited electron into the conduction band (smeared occupations).<br>
         - XCH (Fixed): places the excited electron into the conduction band (fixed occupations).<br>
        <br>
        For XAS calculations of most elements, the FCH treatment is recommended, however in some cases the XCH treatment should be used instead.<br>
        The recommended setting will be shown for each available element.
        Note that only elements for which core-hole pseudopotential sets are available
        will be shown.<br>
        </div>"""
    )
    # I will leave these objects here for now (15/11/23), but since the calculation of molecular
    # systems is not really supported (neither in terms of XAS nor the main App itself) we should
    # not present this option that essentially does nothing.
    # structure_title = ipw.HTML(
    #     """<div style="padding-top: 0px; padding-bottom: 0px">
    #     <h4>Structure</h4></div>"""
    # )
    # structure_help = ipw.HTML(
    #     """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
    #     Below you can indicate if the material should be treated as a molecule
    #     or a crystal.
    #     </div>"""
    # )
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

    def __init__(self, **kwargs):
        self.gipaw_pseudos = pseudo_data_dict["pbe"]["gipaw_pseudos"]
        self.core_hole_pseudos = pseudo_data_dict["pbe"]["core_hole_pseudos"]["1s"]
        self.core_wfc_data_dict = pseudo_data_dict["pbe"]["core_wavefunction_data"]

        self.element_and_ch_treatment = ipw.VBox(layout=ipw.Layout(width="100%"))

        # self.structure_type = ipw.ToggleButtons(
        #     options=[
        #         ("Molecule", "molecule"),
        #         ("Crystal", "crystal"),
        #     ],
        #     value="crystal",
        # )
        self.supercell_min_parameter = ipw.FloatText(
            value=8.0,
            description="The minimum cell length (Å):",
            disabled=False,
            style={"description_width": "initial"},
        )

        self.children = [
            # self.structure_title,
            # self.structure_help,
            # ipw.HBox(
            #     [self.structure_type],
            # ),
            self.element_selection_title,
            self.element_selection_help,
            ipw.HBox([self.element_and_ch_treatment], layout=ipw.Layout(width="95%")),
            self.supercell_title,
            self.supercell_help,
            ipw.HBox(
                [self.supercell_min_parameter],
            ),
        ]

        super().__init__(**kwargs)

    def get_panel_value(self):
        elements_list = []
        core_hole_treatments = {}
        for entry in self.element_and_ch_treatment.children:
            if entry.children[0].value is True:
                element = entry.children[0].description
                ch_treatment = entry.children[1].value
                elements_list.append(element)
                core_hole_treatments[element] = ch_treatment

        pseudo_labels = {}
        core_wfc_data_labels = {}
        for element in elements_list:
            pseudo_labels[element] = {
                "gipaw": self.gipaw_pseudos[element],
                "core_hole": self.core_hole_pseudos[element],
            }
            core_wfc_data_labels[element] = self.core_wfc_data_dict[element]

        parameters = {
            "core_hole_treatments": core_hole_treatments,
            "elements_list": elements_list,
            # "structure_type": self.structure_type.value,
            "pseudo_labels": pseudo_labels,
            "core_wfc_data_labels": core_wfc_data_labels,
            "supercell_min_parameter": self.supercell_min_parameter.value,
        }
        return parameters

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""

        # set selected elements and core-hole treatments
        elements_list = input_dict.get("elements_list", [])
        for entry in self.element_and_ch_treatment.children:
            element = entry.children[0].description
            if element in elements_list:
                entry.children[0].value = True
                entry.children[1].value = input_dict["core_hole_treatments"][element]
            else:
                entry.children[0].value = False
                entry.children[1].value = "full"
        # set supercell min parameter
        self.supercell_min_parameter.value = input_dict.get(
            "supercell_min_parameter", 8.0
        )
        # self.structure_type.value = input_dict.get("structure_type", "crystal")

    @tl.observe("input_structure")
    def _update_structure(self, _=None):
        self._update_element_select_panel()

    def _update_element_select_panel(self):
        if self.input_structure is None:
            return

        starting_treatment_mapping = {"FCH": "full", "XCH": "xch_smear"}
        ch_treatment_options = [
            ("FCH", "full"),
            ("XCH (Smearing)", "xch_smear"),
            ("XCH (Fixed)", "xch_fixed"),
        ]
        ch_pseudos = self.core_hole_pseudos
        structure = self.input_structure
        available_elements = list(ch_pseudos)
        elements_to_select = sorted(
            [
                kind.symbol
                for kind in structure.kinds
                if kind.symbol in available_elements
            ]
        )
        treatment_options = ()

        for element in elements_to_select:
            if element in xch_elements:
                recommended_treatment = "XCH"
            else:
                recommended_treatment = "FCH"

            treatment_options += (
                ipw.HBox(
                    [
                        ipw.Checkbox(
                            description=element,
                            value=False,
                            disabled=False,
                            style={"description_width": "initial"},
                            layout=ipw.Layout(width="7%"),
                        ),
                        ipw.Dropdown(
                            options=ch_treatment_options,
                            value=starting_treatment_mapping[recommended_treatment],
                            disabled=False,
                            layout=ipw.Layout(width="15%"),
                        ),
                        ipw.HTML(
                            f"Recommended treatment: <b>{recommended_treatment}</b> (PBE Core-Hole Pseudopotential)",
                            layout=ipw.Layout(width="78%"),
                        ),
                    ],
                    layout=ipw.Layout(
                        width="100%",
                    ),
                ),
            )

        self.element_and_ch_treatment.children = treatment_options

        # For reference:
        # This is the whole widget:
        # print(f"{self.element_and_ch_treatment}\n")

        # This is the tuple of selected element and core-hole treatment:
        # print(f"{self.element_and_ch_treatment.children[0]}\n")

        # This is the checkbox for the element, giving element name and whether to add it to the elements list
        # print(f"{self.element_and_ch_treatment.children[0].children[0]}\n")
        # print(f"{self.element_and_ch_treatment.children[0].children[0].value}\n")
        # print(f"{self.element_and_ch_treatment.children[0].children[0].description}\n")

        # This is the dropdown for the core-hole treatment option:
        # print(f"{self.element_and_ch_treatment.children[0].children[1]}\n")
        # print(f"{self.element_and_ch_treatment.children[0].children[1].value}\n")

    def reset(self):
        """Reset the panel to its initial state."""
        self.input_structure = None
        # self.structure_type.value = "crystal"
