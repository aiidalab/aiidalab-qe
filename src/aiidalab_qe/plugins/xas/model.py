import os
import tarfile
from importlib import resources
from pathlib import Path

import requests
import traitlets as tl
import yaml

from aiida import orm
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import SettingsModel
from aiidalab_qe.plugins import xas as xas_folder


class XasModel(SettingsModel, HasInputStructure):
    dependencies = [
        "input_structure",
    ]

    # structure_type = tl.Unicode("crystal")
    supercell_min_parameter = tl.Float(8.0)

    kind_names = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
        value_trait=tl.Bool(),  # whether the element is included
    )
    core_hole_treatments = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
        value_trait=tl.Unicode(),  # core hole treatment type
        default_value={},
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        PSEUDO_TOC: dict = yaml.safe_load(
            resources.read_text(
                xas_folder,
                "pseudo_toc.yaml",
            ),
        )  # type: ignore

        self.pseudo_data_dict = PSEUDO_TOC["pseudos"]
        self.xch_elements = PSEUDO_TOC["xas_xch_elements"]
        self.base_url = (
            "https://github.com/PNOGillespie/Core_Level_Spectra_Pseudos/raw/main"
        )
        self.head_path = f"{Path.home()}/.local/lib"
        self.dir_header = "cls_pseudos"
        self.functionals = ["pbe"]
        self.core_wfc_dir = "core_wfc_data"
        self.gipaw_dir = "gipaw_pseudos"
        self.ch_pseudo_dir = "ch_pseudos/star1s"
        self.gipaw_pseudos = self.pseudo_data_dict["pbe"]["gipaw_pseudos"]
        self.core_hole_pseudos = self.pseudo_data_dict["pbe"]["core_hole_pseudos"]["1s"]
        self.core_wfc_data_dict = self.pseudo_data_dict["pbe"]["core_wavefunction_data"]

        self.installed_pseudos = False

    def update(self, specific=""):  # noqa: ARG002
        with self.hold_trait_notifications():
            self._update_pseudos()
            self._update_core_hole_treatment_recommendations()

    def get_model_state(self):
        pseudo_labels = {}
        core_wfc_data_labels = {}
        for element in self.kind_names.keys():
            pseudo_labels[element] = {
                "gipaw": self.gipaw_pseudos[element],
                "core_hole": self.core_hole_pseudos[element],
            }
            core_wfc_data_labels[element] = self.core_wfc_data_dict[element]

        return {
            # "structure_type": self.structure_type,
            "elements_list": list(self.kind_names.keys()),
            "core_hole_treatments": self.core_hole_treatments,
            "pseudo_labels": pseudo_labels,
            "core_wfc_data_labels": core_wfc_data_labels,
            "supercell_min_parameter": self.supercell_min_parameter,
        }

    def set_model_state(self, parameters: dict):
        self.kind_names = {
            kind_name: kind_name in parameters["elements_list"]
            for kind_name in self.kind_names
        }

        self.core_hole_treatments = {
            kind_name: parameters["core_hole_treatments"].get(kind_name, "full")
            for kind_name in self.kind_names
        }

        self.supercell_min_parameter = parameters.get("supercell_min_parameter", 8.0)
        # self.structure_type = parameters.get("structure_type", "crystal")

    def get_kind_names(self):
        return list(self.kind_names)

    def get_recommendation(self, element):
        return "xch_smear" if element in self.xch_elements else "full"

    def reset(self):
        with self.hold_trait_notifications():
            self.supercell_min_parameter = self.traits()[
                "supercell_min_parameter"
            ].default_value
            # self.structure_type = self.traits()["structure_type"].default_value

    def _update_pseudos(self):
        if self.installed_pseudos:
            return

        for func in self.functionals:
            target_dir = f"{self.head_path}/{self.dir_header}/{func}"
            os.makedirs(target_dir, exist_ok=True)
            archive_filename = f"{func}_ch_pseudos.tgz"
            archive_found = any(
                entry == archive_filename for entry in os.listdir(target_dir)
            )
            if not archive_found:
                self._download_extract_pseudo_archive(func)

        # Check all the pseudos/core-wfc data files in the TOC dictionary
        # and load/check all of them before proceeding. Note that this
        # approach relies on there not being multiple instances of nodes
        # with the same label.
        for func in self.functionals:
            gipaw_pseudo_dict = self.pseudo_data_dict[func]["gipaw_pseudos"]
            core_wfc_dict = self.pseudo_data_dict[func]["core_wavefunction_data"]
            core_hole_pseudo_dict = self.pseudo_data_dict[func]["core_hole_pseudos"]
            main_path = f"{self.head_path}/{self.dir_header}/{func}"
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
                self._download_extract_pseudo_archive(func)

            self._load_or_import_nodes_from_filenames(
                in_dict=gipaw_pseudo_dict,
                path=gipaw_dir,
            )
            self._load_or_import_nodes_from_filenames(
                in_dict=core_wfc_dict,
                path=core_wfc_dir,
                core_wfc_data=True,
            )
            self._load_or_import_nodes_from_filenames(
                in_dict=core_hole_pseudo_dict["1s"],
                path=ch_pseudo_dir,
            )

        self.installed_pseudos = True

    def _load_or_import_nodes_from_filenames(self, in_dict, path, core_wfc_data=False):
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

    def _download_extract_pseudo_archive(self, func):
        target_dir = f"{self.head_path}/{self.dir_header}/{func}"
        archive_filename = f"{func}_ch_pseudos.tgz"
        remote_archive_filename = f"{self.base_url}/{func}/{archive_filename}"
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

    def _update_core_hole_treatment_recommendations(self):
        if self.input_structure is None:
            self.kind_names = {}
            self.core_hole_treatments = {}
        else:
            self.kind_names = {
                kind_name: self.kind_names.get(kind_name, False)
                for kind_name in self.input_structure.get_kind_names()
                if kind_name in self.core_hole_pseudos
            }
            self.core_hole_treatments = {
                kind_name: self.get_recommendation(kind_name)
                for kind_name in self.kind_names
            }
