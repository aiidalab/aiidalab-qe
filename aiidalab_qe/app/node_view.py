"""Results view widgets (MOVE TO OTHER MODULE!)

Authors: AiiDAlab team
"""
import base64
import json
import random
import shutil
import typing
from importlib import resources
from pathlib import Path
from tempfile import TemporaryDirectory

import ipywidgets as ipw
import nglview
import numpy as np
import traitlets
from aiida.cmdline.utils.common import get_workchain_report
from aiida.common import LinkType
from aiida.orm import CalcJobNode, Node, ProjectionData, WorkChainNode
from aiidalab_widgets_base import ProcessMonitor, register_viewer_widget
from aiidalab_widgets_base.utils import string_range_to_list
from aiidalab_widgets_base.viewers import StructureDataViewer
from ase import Atoms
from filelock import FileLock, Timeout
from IPython.display import HTML, clear_output, display
from jinja2 import Environment
from monty.json import jsanitize
from traitlets import Instance, Int, List, Unicode, Union, default, observe, validate
from widget_bandsplot import BandsPlotWidget

from aiidalab_qe.app import static
from aiidalab_qe.app.report import generate_report_html, get_hubbard_occupations_list


class MinimalStructureViewer(ipw.VBox):
    structure = Union([Instance(Atoms), Instance(Node)], allow_none=True)
    _displayed_structure = Instance(Atoms, allow_none=True, read_only=True)

    background = Unicode()
    supercell = List(Int())

    def __init__(self, structure, *args, **kwargs):
        self._viewer = nglview.NGLWidget()
        self._viewer.camera = "orthographic"
        self._viewer.stage.set_parameters(mouse_preset="pymol")
        ipw.link((self, "background"), (self._viewer, "background"))

        self.structure = structure

        super().__init__(
            children=[
                self._viewer,
            ],
            *args,
            **kwargs,
        )

    @default("background")
    def _default_background(self):
        return "#FFFFFF"

    @default("supercell")
    def _default_supercell(self):
        return [1, 1, 1]

    @validate("structure")
    def _valid_structure(self, change):  # pylint: disable=no-self-use
        """Update structure."""
        structure = change["value"]

        if structure is None:
            return None  # if no structure provided, the rest of the code can be skipped

        if isinstance(structure, Atoms):
            return structure
        if isinstance(structure, Node):
            return structure.get_ase()
        raise ValueError(
            "Unsupported type {}, structure must be one of the following types: "
            "ASE Atoms object, AiiDA CifData or StructureData."
        )

    @observe("structure")
    def _update_displayed_structure(self, change):
        """Update displayed_structure trait after the structure trait has been modified."""
        # Remove the current structure(s) from the viewer.
        if change["new"] is not None:
            self.set_trait("_displayed_structure", change["new"].repeat(self.supercell))
        else:
            self.set_trait("_displayed_structure", None)

    @observe("_displayed_structure")
    def _update_structure_viewer(self, change):
        """Update the view if displayed_structure trait was modified."""
        with self.hold_trait_notifications():
            for (
                comp_id
            ) in self._viewer._ngl_component_ids:  # pylint: disable=protected-access
                self._viewer.remove_component(comp_id)
            self.selection = list()
            if change["new"] is not None:
                self._viewer.add_component(nglview.ASEStructure(change["new"]))
                self._viewer.clear()
                self._viewer.stage.set_parameters(clipDist=0)
                self._viewer.add_representation("unitcell", diffuse="#df0587")
                self._viewer.add_representation("ball+stick", aspectRatio=3.5)


def export_bands_data(work_chain_node, fermi_energy=None):
    if "band_structure" in work_chain_node.outputs:
        data = json.loads(
            work_chain_node.outputs.band_structure._exportcontent(
                "json", comments=False
            )[0]
        )
        # The fermi energy from band calculation is not robust.
        data["fermi_level"] = (
            fermi_energy or work_chain_node.outputs.band_parameters["fermi_energy"]
        )
        return [
            jsanitize(data),
        ]
    else:
        return None


def cmap(label: str) -> str:
    """Return RGB string of color for given pseudo info
    Hardcoded at the momment.
    """
    # if a unknow type generate random color based on ascii sum
    ascn = sum([ord(c) for c in label])
    random.seed(ascn)

    return "#%06x" % random.randint(0, 0xFFFFFF)


def _projections_curated(
    projections: ProjectionData,
    group_dos_by="atom",
    spin_type="none",
    line_style="solid",
):
    """Collect the data from ProjectionData and parse it as dos list which can be
    understand by bandsplot widget. `group_dos_by` is for which tag to be grouped, by atom or by orbital name.
    The spin_type is used to invert all the y values of pdos to be shown as spin down pdos and to set label.
    """
    _pdos = {}

    for orbital, pdos, energy in projections.get_pdos():
        orbital_data = orbital.get_orbital_dict()
        kind_name = orbital_data["kind_name"]
        atom_position = [round(i, 2) for i in orbital_data["position"]]
        orbital_name = orbital.get_name_from_quantum_numbers(
            orbital_data["angular_momentum"], orbital_data["magnetic_number"]
        ).lower()

        if group_dos_by == "atom":
            dos_group_name = atom_position
        elif group_dos_by == "angular":
            # by orbital label
            dos_group_name = orbital_name[0]
        elif group_dos_by == "angular_and_magnetic":
            # by orbital label
            dos_group_name = orbital_name
        else:
            raise Exception(f"Unknow dos type: {group_dos_by}!")

        key = f"{kind_name}-{dos_group_name}"
        if key in _pdos:
            _pdos[key][1] += pdos
        else:
            _pdos[key] = [energy, pdos]

    dos = []
    for label, (energy, pdos) in _pdos.items():
        if spin_type == "down":
            # invert y-axis
            pdos = -pdos
            label = f"{label} (↓)"

        if spin_type == "up":
            label = f"{label} (↑)"

        orbital_pdos = {
            "label": label,
            "x": energy.tolist(),
            "y": pdos.tolist(),
            "borderColor": cmap(label),
            "lineStyle": line_style,
        }
        dos.append(orbital_pdos)

    return dos


def export_pdos_data(work_chain_node, group_dos_by="atom"):
    if "dos" in work_chain_node.outputs:
        _, energy_dos, _ = work_chain_node.outputs.dos.get_x()
        tdos_values = {f"{n}": v for n, v, _ in work_chain_node.outputs.dos.get_y()}

        dos = []

        if "projections" in work_chain_node.outputs:
            # The total dos parsed
            tdos = {
                "label": "Total DOS",
                "x": energy_dos.tolist(),
                "y": tdos_values.get("dos").tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            dos.append(tdos)

            dos += _projections_curated(
                work_chain_node.outputs.projections,
                group_dos_by=group_dos_by,
                spin_type="none",
            )

        else:
            # The total dos parsed
            tdos_up = {
                "label": "Total DOS (↑)",
                "x": energy_dos.tolist(),
                "y": tdos_values.get("dos_spin_up").tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            tdos_down = {
                "label": "Total DOS (↓)",
                "x": energy_dos.tolist(),
                "y": (-tdos_values.get("dos_spin_down")).tolist(),  # minus
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "dash",
            }
            dos += [tdos_up, tdos_down]

            # spin-up (↑)
            dos += _projections_curated(
                work_chain_node.outputs.projections_up,
                group_dos_by=group_dos_by,
                spin_type="up",
            )

            # spin-dn (↓)
            dos += _projections_curated(
                work_chain_node.outputs.projections_down,
                group_dos_by=group_dos_by,
                spin_type="down",
                line_style="dash",
            )

        data_dict = {
            "fermi_energy": work_chain_node.outputs.nscf_parameters["fermi_energy"],
            "dos": dos,
        }

        return json.loads(json.dumps(data_dict))

    else:
        return None


def export_data(work_chain_node, group_dos_by="atom"):
    dos = export_pdos_data(work_chain_node, group_dos_by=group_dos_by)
    fermi_energy = dos["fermi_energy"] if dos else None

    bands = export_bands_data(work_chain_node, fermi_energy)

    return dict(
        bands=bands,
        dos=dos,
    )


class VBoxWithCaption(ipw.VBox):
    def __init__(self, caption, body, *args, **kwargs):
        super().__init__(children=[ipw.HTML(caption), body], *args, **kwargs)


class SummaryView(ipw.VBox):
    def __init__(self, wc_node, **kwargs):
        report_html = generate_report_html(wc_node)

        self.summary_view = ipw.HTML(report_html)
        super().__init__(
            children=[self.summary_view],
            **kwargs,
        )


class WorkChainOutputs(ipw.VBox):
    _busy = traitlets.Bool(read_only=True)

    def __init__(self, node, export_dir=None, **kwargs):
        self.export_dir = Path.cwd().joinpath("exports")

        if node.process_label != "QeAppWorkChain":
            raise KeyError(str(node.node_type))

        self.node = node

        self._create_archive_indicator = ipw.HTML(
            """<button disabled>
                <i class="fa fa-spinner fa-spin" aria-hidden="true"></i>
                    Creating archive...
            </button>"""
        )
        self._download_archive_button = ipw.Button(
            description="Download output",
            icon="download",
        )
        self._download_archive_button.on_click(self._download_archive)
        self._download_button_container = ipw.Box([self._download_archive_button])

        if node.exit_status != 0:
            title = ipw.HTML(
                f"<h4>Workflow failed with exit status [{ node.exit_status }]</h4>"
            )
            final_calcjob = self._get_final_calcjob(node)
            env = Environment()
            template = resources.read_text(static, "workflow_failure.jinja")
            style = resources.read_text(static, "style.css")
            output = ipw.HTML(
                env.from_string(template).render(
                    style=style,
                    process_report=get_workchain_report(node, "REPORT"),
                    calcjob_exit_message=final_calcjob.exit_message,
                )
            )
        else:
            title = ipw.HTML("<h4>Workflow completed successfully!</h4>")
            output = ipw.HTML()

        super().__init__(
            children=[
                ipw.HBox(
                    children=[title, self._download_button_container],
                    layout=ipw.Layout(justify_content="space-between", margin="10px"),
                ),
                output,
            ],
            **kwargs,
        )

    @traitlets.default("_busy")
    def _default_busy(self):
        return False

    @traitlets.observe("_busy")
    def _observe_busy(self, change):
        self._download_button_container.children = [
            self._create_archive_indicator
            if change["new"]
            else self._download_archive_button
        ]

    def _download_archive(self, _):
        fn_archive = self.export_dir.joinpath(str(self.node.uuid)).with_suffix(".zip")
        fn_lockfile = fn_archive.with_suffix(".lock")

        try:
            self.set_trait("_busy", True)
            # Create exports archive directory.
            fn_archive.parent.mkdir(parents=True, exist_ok=True)
            # Try to obtain lock for creating archive...
            with FileLock(fn_lockfile, timeout=0):
                # Check whether archive file already exists.
                if not fn_archive.is_file():
                    # Create archive file.
                    with TemporaryDirectory() as tmpdir:
                        self._prepare_calcjob_io(self.node, Path(tmpdir))
                        shutil.make_archive(fn_archive.with_suffix(""), "zip", tmpdir)
            Path(fn_lockfile).unlink()  # Delete lock file.
        except Timeout:
            # Failed to obtain lock, presuming some other process is working on it.
            with FileLock(fn_lockfile, timeout=20):
                assert fn_archive.is_file()
        finally:
            self.set_trait("_busy", False)

        id = f"dl_{self.node.uuid}"

        display(
            HTML(
                f"""
        <html>
        <body>
        <a
            id="{id}"
            href="{fn_archive.relative_to(Path.cwd())}"
            download="{fn_archive.stem}"
        ></a>
        <script>
        (function download() {{document.getElementById("{id}").click();
        }})()
        </script>

        </body>
        </html>
        """
            )
        )

    @classmethod
    def _prepare_calcjob_io(cls, node: WorkChainNode, root_folder: Path):
        """Prepare the calculation job input and output files.

        :param node: QeAppWorkChain node.
        """
        counter = 1

        for link1 in node.get_outgoing(link_type=LinkType.CALL_WORK):
            wc_node = link1.node
            for link2 in wc_node.get_outgoing(link_type=LinkType.CALL_WORK):
                base_node = link2.node
                base_label = (
                    f"iter{link2.link_label[-1]}"
                    if link1.link_label == "relax"
                    else link2.link_label
                )
                for link3 in base_node.get_outgoing(
                    link_type=LinkType.CALL_CALC, link_label_filter="iteration_%"
                ):
                    counter_str = f"0{counter}" if counter < 10 else str(counter)
                    pw_label = f"pw{link3.link_label[-1]}"

                    fdname = f"{counter_str}-{link1.link_label}-{base_label}-{pw_label}"

                    folder_path = root_folder / fdname

                    cls._write_calcjob_io(link3.node, folder_path)

                    counter += 1

    @staticmethod
    def _get_final_calcjob(node: WorkChainNode) -> typing.Union[None, CalcJobNode]:
        """Get the final calculation job node called by a work chain node.

        :param node: Work chain node.
        """
        try:
            final_calcjob = [
                process
                for process in node.called_descendants
                if isinstance(process, CalcJobNode) and process.is_finished
            ][-1]
        except IndexError:
            final_calcjob = None

        return final_calcjob

    @staticmethod
    def _write_calcjob_io(calcjob: CalcJobNode, folder: Path) -> None:
        """Write the ``calcjob`` in and output files to ``folder``.

        :param calcjob: calculation job node for which to write the IO files.
        :param folder: folder to which to write the IO files.
        """
        folder.mkdir(exist_ok=True)
        input_filepath = folder / "aiida.in"

        with calcjob.open(calcjob.get_option("input_filename"), "r") as ihandle:
            with input_filepath.open("w") as ohandle:
                ohandle.write(ihandle.read())

        pseudo_folder = folder / "pseudo"
        pseudo_folder.mkdir(exist_ok=True)

        for _, pseudo in calcjob.inputs.pseudos.items():
            pseudo_path = pseudo_folder / pseudo.filename

            with pseudo_path.open("w") as handle:
                handle.write(pseudo.get_content())

        retrieved = calcjob.outputs.retrieved

        for filename in retrieved.list_object_names():
            out_filepath = folder / filename
            with out_filepath.open("w") as handle:
                handle.write(retrieved.get_object_content(filename))


@register_viewer_widget("process.workflow.workchain.WorkChainNode.")
class WorkChainViewer(ipw.VBox):
    _results_shown = traitlets.Set()

    def __init__(self, node, **kwargs):
        if node.process_label != "QeAppWorkChain":
            super().__init__()
            return

        self.node = node

        self.title = ipw.HTML(
            f"""
            <hr style="height:2px;background-color:#2097F3;">
            <h4>QE App Workflow (pk: {self.node.pk}) &mdash;
                {self.node.inputs.structure.get_formula()}
            </h4>
            """
        )
        self.workflows_summary = SummaryView(self.node)

        self.summary_tab = ipw.VBox(children=[self.workflows_summary])
        self.structure_tab = ipw.VBox(
            [ipw.Label("Structure not available.")],
            layout=ipw.Layout(min_height="380px"),
        )
        self.bands_tab = ipw.VBox(
            [ipw.Label("Electronic Structure not available.")],
            layout=ipw.Layout(min_height="380px"),
        )

        # Check if DFT+U Calculation was conducted
        dict_parameters = self.node.base.extras.get("builder_parameters", {})
        self.hubbard_cond = dict_parameters.get("hubbard")

        if self.hubbard_cond == "Yes":
            self.hubbard_tab = ipw.VBox(
                [ipw.Label("Hubbard occupations not available.")],
                layout=ipw.Layout(min_height="380px"),
            )
            self.result_tabs = ipw.Tab(
                children=[
                    self.summary_tab,
                    self.structure_tab,
                    self.bands_tab,
                    self.hubbard_tab,
                ]
            )

            self.result_tabs.set_title(0, "Workflow Summary")
            self.result_tabs.set_title(1, "Final Geometry (n/a)")
            self.result_tabs.set_title(2, "Electronic Structure (n/a)")
            self.result_tabs.set_title(3, "DFT+U (n/a)")
        else:
            self.result_tabs = ipw.Tab(
                children=[self.summary_tab, self.structure_tab, self.bands_tab]
            )
            self.result_tabs.set_title(0, "Workflow Summary")
            self.result_tabs.set_title(1, "Final Geometry (n/a)")
            self.result_tabs.set_title(2, "Electronic Structure (n/a)")

        # An ugly fix to the structure appearance problem
        # https://github.com/aiidalab/aiidalab-qe/issues/69
        def on_selected_index_change(change):
            index = change["new"]
            # Accessing the viewer only if the corresponding tab is present.
            if self.result_tabs._titles[str(index)] == "Final Geometry":
                self._structure_view._viewer.handle_resize()

                def toggle_camera():
                    """Toggle camera between perspective and orthographic."""
                    self._structure_view._viewer.camera = (
                        "perspective"
                        if self._structure_view._viewer.camera == "orthographic"
                        else "orthographic"
                    )

                toggle_camera()
                toggle_camera()

        self.result_tabs.observe(on_selected_index_change, "selected_index")
        self._update_view()

        super().__init__(
            children=[self.title, self.result_tabs],
            **kwargs,
        )
        self._process_monitor = ProcessMonitor(
            process=self.node,
            callbacks=[
                self._update_view,
            ],
        )

    def _update_view(self):
        with self.hold_trait_notifications():
            if self.node.is_finished:
                self._show_workflow_output()
                if self.node.is_finished_ok and self.hubbard_cond == "Yes":
                    self._show_hubbard_occupations()
            if (
                "structure" not in self._results_shown
                and "structure" in self.node.outputs
            ):
                self._show_structure()
                self._results_shown.add("structure")

            if "electronic_structure" not in self._results_shown and (
                "band_structure" in self.node.outputs or "dos" in self.node.outputs
            ):
                self._show_electronic_dosoptions()
                # self._show_electronic_structure()
                self._results_shown.add("electronic_structure")

    def _show_structure(self):
        self._structure_view = StructureDataViewer(
            structure=self.node.outputs.structure
        )
        self.result_tabs.children[1].children = [self._structure_view]
        self.result_tabs.set_title(1, "Final Geometry")

    def _show_electronic_dosoptions(self):
        self.dos_options = BandDosPlotsWidget(self.node)
        self.result_tabs.children[2].children = [self.dos_options]
        self.result_tabs.set_title(2, "Electronic Structure")

    def _show_electronic_structure(self):
        group_dos_by = ipw.ToggleButtons(
            options=[
                ("Atom", "atom"),
                ("Orbital", "angular"),
            ],
            value="atom",
        )
        settings = ipw.VBox(
            children=[
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "DOS grouped by:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="120px"
                            ),
                        ),
                        group_dos_by,
                    ]
                ),
            ],
            layout={"margin": "0 0 30px 30px"},
        )
        #
        data = export_data(self.node, group_dos_by=group_dos_by.value)
        bands_data = data.get("bands", None)
        dos_data = data.get("dos", None)
        _bands_plot_view = BandsPlotWidget(
            bands=bands_data,
            dos=dos_data,
        )

        def response(change):
            data = export_data(self.node, group_dos_by=group_dos_by.value)
            bands_data = data.get("bands", None)
            dos_data = data.get("dos", None)
            _bands_plot_view = BandsPlotWidget(
                bands=bands_data,
                dos=dos_data,
            )
            self.result_tabs.children[2].children = [
                settings,
                _bands_plot_view,
            ]

        group_dos_by.observe(response, names="value")
        # update the electronic structure tab
        self.result_tabs.children[2].children = [
            settings,
            _bands_plot_view,
        ]
        self.result_tabs.set_title(2, "Electronic Structure")

    def _show_workflow_output(self):
        self.workflows_output = WorkChainOutputs(self.node)

        self.result_tabs.children[0].children = [
            self.workflows_summary,
            self.workflows_output,
        ]

    def _show_hubbard_occupations(self):
        self.hubbard_summary = HubbardView(self.node)
        self.result_tabs.children[3].children = [self.hubbard_summary]
        self.result_tabs.set_title(3, "DFT+U")


class HubbardView(ipw.VBox):
    def __init__(self, node, **kwargs):
        self.node = node
        self.data_occupations = get_hubbard_occupations_list(self.node)
        self.hubbard_keys = self._get_hubbard_keys()
        self.hubbard_options = self._get_options_dict()
        self.spin_type = self._get_spin_type()
        self.select_kind = ipw.Dropdown(options=self.hubbard_keys, description="Kind:")
        self.select_index = ipw.Dropdown(
            options=self.hubbard_options[self.select_kind.value], description="Atom:"
        )
        self.env = Environment()
        self.template = resources.read_text(static, "hubbard_occupation.jinja")
        self.style = resources.read_text(static, "style.css")
        self.hubbard_occupation = ipw.HTML(
            self.env.from_string(self.template).render(
                style=self.style, **self._select_dict()
            )
        )

        super().__init__(
            children=[
                ipw.HBox([self.select_kind, self.select_index]),
                self.hubbard_occupation,
            ],
        ),
        self.select_kind.observe(self._update_index, names="value")
        self.select_index.observe(self._update_hubbard_occupation, names="value")

    def _get_spin_type(self):
        dict_parameters = self.node.base.extras.get("builder_parameters")
        input_spin_type = dict_parameters.get("spin_type")
        if input_spin_type == "none":
            spin_type = False
        else:
            spin_type = True
        return spin_type

    def _update_index(self, change):
        self.select_index.options = self.hubbard_options[change["new"]]

    def _select_dict(self):
        data = next(
            (
                d
                for d in self.data_occupations
                if d.get("atom_index") == self.select_index.value
            ),
            {},
        )
        data["spin_type"] = self.spin_type
        return data

    def _update_hubbard_occupation(self, change):
        data = next(
            (d for d in self.data_occupations if d.get("atom_index") == change["new"]),
            None,
        )
        data["spin_type"] = self.spin_type
        self.hubbard_occupation.value = self.env.from_string(self.template).render(
            style=self.style, **data
        )

    def _get_hubbard_keys(self):
        dict_parameters = self.node.base.extras.get("builder_parameters")
        hubbard_dict = dict_parameters.get("hubbard_dict", None)
        return list(hubbard_dict.keys())

    def _get_options_dict(self):
        list_indixes = []
        sites = self.node.inputs.structure.sites
        for kind in self.hubbard_keys:
            list_temp = [
                str(index + 1)
                for index, site in enumerate(sites)
                if site.kind_name == kind
            ]
            list_indixes.append(list_temp)
        options_dict = dict(zip(self.hubbard_keys, list_indixes))
        return options_dict


class BandDosPlotsWidget(ipw.VBox):
    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Select the style of plotting the projected density of states.
        </div>"""
    )

    def __init__(self, node, **kwargs):
        self.node = node

        self.dos_atoms_group = ipw.Dropdown(
            options=[("Kinds", "kinds"), ("Atoms", "atoms")],
            value="kinds",
        )
        self.dos_plot_group = ipw.Dropdown(
            options=[("Total", "total"), ("Orbital", "orbital")],
            value="total",
        )
        self.selected_atoms = ipw.Text(
            description="Select atoms:",
            value="",
            style={"description_width": "initial"},
        )
        self.wrong_syntax = ipw.HTML(
            value="""<i class="fa fa-times" style="color:red;font-size:2em;" ></i> wrong syntax""",
            layout={"visibility": "hidden"},
        )
        self.update_plot_button = ipw.Button(
            description="Update Plot",
            icon="pencil",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )
        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto", visibility="hidden"),
        )
        self.export_dir = Path.cwd().joinpath("exports")
        self.dos_data = self._get_dos_data()
        self.fermi_energy = self._get_fermi_energy()
        self.bands_data = self._get_bands_data()
        self.band_labels = self.bands_labeling(self.bands_data)
        self.bands_xaxis = self._band_xaxis()
        self.bands_yaxis = self._band_yaxis()
        self.dos_xaxis = self._dos_xaxis()
        self.dos_yaxis = self._dos_yaxis()
        self.bandsplot_widget = self._bandsplot_widget()

        if self.bands_data and not self.dos_data:
            self.update_plot_button.disabled = True
        self.bands_widget = ipw.Output()

        def download_data(_=None):
            file_name_bands = "bands_data.json"
            file_name_dos = "dos_data.json"
            if self.bands_data:
                json_str = json.dumps(self.bands_data)
                b64_str = base64.b64encode(json_str.encode()).decode()
                self._download(payload=b64_str, filename=file_name_bands)
            if self.dos_data:
                json_str = json.dumps(self.dos_data)
                b64_str = base64.b64encode(json_str.encode()).decode()
                self._download(payload=b64_str, filename=file_name_dos)

        self.download_button.on_click(download_data)
        self._initial_view()
        self.update_plot_button.on_click(self._update_plot)

        super().__init__(
            children=[
                self.description,
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Group :",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="80px"
                            ),
                        ),
                        self.dos_atoms_group,
                    ]
                ),
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Plot :",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="80px"
                            ),
                        ),
                        self.dos_plot_group,
                    ]
                ),
                ipw.HBox([self.selected_atoms, self.wrong_syntax]),
                ipw.HBox(
                    children=[
                        self.update_plot_button,
                        self.download_button,
                    ]
                ),
                self.bands_widget,
            ]
        )

    @staticmethod
    def _download(payload, filename):
        """Download payload as a file named as filename."""
        from IPython.display import Javascript

        javas = Javascript(
            """
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """.format(
                payload=payload, filename=filename
            )
        )
        display(javas)

    def _get_dos_data(self):
        if "pdos" in self.node.inputs.properties:
            expanded_selection, syntax_ok = string_range_to_list(
                self.selected_atoms.value, shift=-1
            )
            dos = get_pdos_data(
                self.node,
                group_tag=self.dos_atoms_group.value,
                plot_tag=self.dos_plot_group.value,
                selected_atoms=expanded_selection,
            )
            return dos
        else:
            return None

    def _get_fermi_energy(self):
        fermi_energy = self.dos_data["fermi_energy"] if self.dos_data else 0.0
        return fermi_energy

    def _get_bands_data(self):
        if "bands" in self.node.inputs.properties:
            bands = export_bands_data(self.node, self.fermi_energy)
            return bands
        else:
            return None

    def _initial_view(self):
        with self.bands_widget:
            clear_output(wait=True)
            # self.bandsplot_widget.show() #Fix plotly not showing
            display(self.bandsplot_widget)
            self.download_button.layout.visibility = "visible"

    def _update_plot(self, _=None):
        with self.bands_widget:
            expanded_selection, syntax_ok = string_range_to_list(
                self.selected_atoms.value, shift=-1
            )
            if not syntax_ok:
                self.wrong_syntax.layout.visibility = "visible"
                clear_output(wait=True)
            else:
                self.dos_data = self._get_dos_data()
                self.bandsplot_widget = self._bandsplot_widget()
                clear_output(wait=True)
                # self.bandsplot_widget.show() #Fix plotly not showing
                display(self.bandsplot_widget)

    def _bandsplot_widget(self):
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        if self.bands_data and not self.dos_data:
            fig = go.Figure()
            paths = self.bands_data[0].get("paths")

            for band in paths:
                for bands in band["values"]:
                    bands_np = np.array(bands)
                    fig.add_trace(
                        go.Scatter(
                            x=band["x"],
                            y=bands_np - self.fermi_energy,  # substract Fermi Energy
                            mode="lines",
                            line=dict(color="#111111", shape="spline", smoothing=1.3),
                            showlegend=False,
                        )
                    )
            for i in self.band_labels[1]:
                fig.add_vline(x=i, line=dict(color="#111111", width=1))
            if self.fermi_energy != 0:
                fig.add_hline(y=0, line=dict(color="#111111", width=1, dash="dot"))
            fig.update_layout(
                height=600,
                width=950,
                plot_bgcolor="white",
                xaxis=self.bands_xaxis,
                yaxis=self.bands_yaxis,
            )

        elif self.dos_data and not self.bands_data:
            fig = go.Figure()
            for trace in self.dos_data["dos"]:
                if trace["label"] == "Total DOS":
                    my_fill = "tozeroy"
                else:
                    my_fill = None
                dos_np = np.array(trace["x"])
                fig.add_trace(
                    go.Scatter(
                        x=dos_np - self.fermi_energy,
                        y=trace["y"],
                        fill=my_fill,
                        name=trace["label"],
                        line=dict(
                            color=trace["borderColor"], shape="spline", smoothing=1.0
                        ),
                    )
                )
            fig.add_vline(x=0, line=dict(color="#111111", width=1, dash="dot"))
            fig.update_layout(
                height=600,
                width=850,
                plot_bgcolor="white",
                xaxis=self.dos_xaxis,
                yaxis=self.dos_yaxis,
            )

        elif self.dos_data and self.bands_data:
            fig = make_subplots(
                rows=1,
                cols=2,
                shared_yaxes=True,
                column_widths=[0.7, 0.3],
                horizontal_spacing=0.02,
            )
            paths = self.bands_data[0].get("paths")
            for band in paths:
                for bands in band["values"]:
                    bands_np = np.array(bands)
                    fig.add_trace(
                        go.Scatter(
                            x=band["x"],
                            y=bands_np - self.fermi_energy,  # substract Fermi Energy
                            mode="lines",
                            line=dict(color="#111111", shape="spline", smoothing=1.3),
                            showlegend=False,
                        ),
                        row=1,
                        col=1,
                    )

            for trace in self.dos_data["dos"]:
                if trace["label"] == "Total DOS":
                    my_fill = "tozerox"
                else:
                    my_fill = None

                dos_np = np.array(trace["x"])
                fig.add_trace(
                    go.Scatter(
                        x=trace["y"],
                        y=dos_np - self.fermi_energy,
                        fill=my_fill,
                        name=trace["label"],
                        line=dict(
                            color=trace["borderColor"], shape="spline", smoothing=1.3
                        ),
                    ),
                    row=1,
                    col=2,
                )
            for i in self.band_labels[1]:
                fig.add_vline(
                    x=i,
                    line=dict(color="#111111", width=1),
                    row=1,
                    col=1,
                )
            fig.update_xaxes(
                patch=self.bands_xaxis,
                row=1,
                col=1,
            )
            fig.update_yaxes(
                patch=self.bands_yaxis,
                row=1,
                col=1,
            )
            fig.update_xaxes(
                patch=self.dos_xaxis,
                row=1,
                col=2,
            )
            fig.update_yaxes(patch=self.dos_yaxis, row=1, col=2, showticklabels=False)

            fig.add_hline(
                y=0, line=dict(color="#111111", width=1, dash="dot"), row=1, col=1
            )
            fig.add_hline(
                y=0, line=dict(color="#111111", width=1, dash="dot"), row=1, col=2
            )
            fig.update_layout(
                height=600,
                width=850,
                plot_bgcolor="white",
                legend=dict(
                    # yanchor="top",
                    # y=0.99,
                    xanchor="left",
                    x=1.04,
                ),
            )

        else:
            fig = None

        return go.FigureWidget(fig)

    def _band_xaxis(self):
        import plotly.graph_objects as go

        if self.bands_data:
            paths = self.bands_data[0].get("paths")
            # labels = bands_labeling(self.bands_data)
            # labels, labels_values = bands_labeling(self.bands_data)
            slider_bands = go.layout.xaxis.Rangeslider(
                thickness=0.08,
                range=[0, paths[-1]["x"][-1]],
            )

            bandxaxis = go.layout.XAxis(
                title="k-points",
                range=[0, paths[-1]["x"][-1]],
                showgrid=True,
                showline=True,
                tickmode="array",
                rangeslider=slider_bands,
                fixedrange=False,
                tickvals=self.band_labels[1],
                ticktext=self.band_labels[0],
                showticklabels=True,
                linecolor="#111111",
                mirror=True,
                linewidth=2,
                type="linear",
            )

        else:
            bandxaxis = None
        return bandxaxis

    def _band_yaxis(self):
        import plotly.graph_objects as go

        if self.bands_data:
            bandyaxis = go.layout.YAxis(
                # title="Electronic Bands(eV)",
                title=dict(text="Electronic Bands (eV)", standoff=1),
                side="left",
                showgrid=True,
                showline=True,
                zeroline=True,
                fixedrange=False,
                automargin=True,
                mirror="ticks",
                ticks="inside",
                linewidth=2,
                linecolor="#111111",
                tickwidth=2,
                zerolinewidth=2,
            )

        else:
            bandyaxis = None
        return bandyaxis

    def _dos_xaxis(self):
        import plotly.graph_objects as go

        if self.dos_data:
            if self.bands_data:
                dosxaxis = go.layout.XAxis(
                    title="Density of states",
                    side="bottom",
                    showgrid=True,
                    showline=True,
                    linecolor="#111111",
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                )

            else:
                dosxaxis = go.layout.XAxis(
                    title="Density of states (eV)",
                    showgrid=True,
                    showline=True,
                    linecolor="#111111",
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                )

        else:
            dosxaxis = None
        return dosxaxis

    def _dos_yaxis(self):
        import plotly.graph_objects as go

        if self.dos_data:
            if self.bands_data:
                dosyaxis = go.layout.YAxis(
                    # title="Density of states (eV)", #FigureWidget modifies the title position in the meantime lets put it
                    showgrid=True,
                    showline=True,
                    side="right",
                    # position=0.0,
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                    linecolor="#111111",
                    zerolinewidth=2,
                )

            else:
                dosyaxis = go.layout.YAxis(
                    # title="Density of states (eV)",
                    showgrid=True,
                    showline=True,
                    side="left",
                    mirror="ticks",
                    ticks="inside",
                    linewidth=2,
                    tickwidth=2,
                    linecolor="#111111",
                    zerolinewidth=2,
                )

        else:
            dosyaxis = None
        return dosyaxis

    def bands_labeling(self, bands):
        """Function to return two list containing the labels and values (kpoint) for plotting
        params: bands: dictionary from export_bands_data function
        output: label (list of str), label_values (list of float)
        """
        if bands:
            my_paths = bands[0].get("paths")
            my_labels = []
            for path in my_paths:  # Remove duplicates
                label_a = [path["from"], path["x"][0]]
                label_b = [path["to"], path["x"][-1]]
                if label_a not in my_labels:
                    my_labels.append(label_a)
                if label_b not in my_labels:
                    my_labels.append(label_b)

            my_clean_labels = []  # Format
            for i in my_labels:
                if my_clean_labels:
                    if i not in my_clean_labels:
                        if my_clean_labels[-1][-1] == i[1]:
                            my_clean_labels[-1][0] = my_clean_labels[-1][0] + "|" + i[0]
                        else:
                            my_clean_labels.append(i)
                else:
                    my_clean_labels.append(i)

            labels = [label[0] for label in my_clean_labels]
            labels_values = [label[1] for label in my_clean_labels]
            return [labels, labels_values]
        else:
            return None


def get_pdos_data(work_chain_node, group_tag, plot_tag, selected_atoms):
    if "dos" in work_chain_node.outputs:
        _, energy_dos, _ = work_chain_node.outputs.dos.get_x()
        tdos_values = {f"{n}": v for n, v, _ in work_chain_node.outputs.dos.get_y()}

        dos = []

        if "projections" in work_chain_node.outputs:
            # The total dos parsed
            tdos = {
                "label": "Total DOS",
                "x": energy_dos.tolist(),
                "y": tdos_values.get("dos").tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            dos.append(tdos)

            dos += _projections_curated_options(
                work_chain_node.outputs.projections,
                spin_type="none",
                group_tag=group_tag,
                plot_tag=plot_tag,
                selected_atoms=selected_atoms,
            )

        else:
            # The total dos parsed
            tdos_up = {
                "label": "Total DOS (↑)",
                "x": energy_dos.tolist(),
                "y": tdos_values.get("dos_spin_up").tolist(),
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "solid",
            }
            tdos_down = {
                "label": "Total DOS (↓)",
                "x": energy_dos.tolist(),
                "y": (-tdos_values.get("dos_spin_down")).tolist(),  # minus
                "borderColor": "#8A8A8A",  # dark gray
                "backgroundColor": "#999999",  # light gray
                "backgroundAlpha": "40%",
                "lineStyle": "dash",
            }
            dos += [tdos_up, tdos_down]

            # spin-up (↑)
            dos += _projections_curated_options(
                work_chain_node.outputs.projections_up,
                spin_type="up",
                group_tag=group_tag,
                plot_tag=plot_tag,
                selected_atoms=selected_atoms,
            )

            # spin-dn (↓)
            dos += _projections_curated_options(
                work_chain_node.outputs.projections_down,
                spin_type="down",
                line_style="dash",
                group_tag=group_tag,
                plot_tag=plot_tag,
                selected_atoms=selected_atoms,
            )

        data_dict = {
            "fermi_energy": work_chain_node.outputs.nscf_parameters["fermi_energy"],
            "dos": dos,
        }

        return json.loads(json.dumps(data_dict))

    else:
        return None


def _projections_curated_options(
    projections: ProjectionData,
    group_tag,
    plot_tag,
    selected_atoms,
    spin_type="none",
    line_style="solid",
):
    """Collect the data from ProjectionData and parse it as dos list which can be
    understand by bandsplot widget. `curated_tag` is for which tag to be grouped, by atom or by orbital name.
    The spin_type is used to invert all the y values of pdos to be shown as spin down pdos and to set label.
    """
    _pdos = {}
    list_positions = []
    # Setting
    dict_html = {
        "pz": "p<sub>z</sub>",
        "px": "p<sub>x</sub>",
        "py": "p<sub>y</sub>",
        "dz2": "d<sub>z<sup>2</sup></sub>",
        "dxy": "d<sub>xy</sub>",
        "dxz": "d<sub>xz</sub>",
        "dyz": "d<sub>yz</sub>",
        "dx2-y2": "d<sub>x<sup>2</sup>-y<sup>2</sup></sub>",
        "fz3": "f<sub>z<sup>3</sup></sub>",
        "fxz2": "f<sub>xz<sup>2</sup></sub>",
        "fyz2": "f<sub>yz<sup>2</sup></sub>",
        "fxyz": "f<sub>xzy</sub>",
        "fx(x2-3y2)": "f<sub>x(x<sup>2</sup>-3y<sup>2</sup>)</sub>",
        "fy(3x2-y2)": "f<sub>y(3x<sup>2</sup>-y<sup>2</sup>)</sub>",
        "fy(x2-z2)": "f<sub>y(x<sup>2</sup>-z<sup>2</sup>)</sub>",
        0.5: "<sup>+1</sup>/<sub>2</sub>",
        -0.5: "<sup>-1</sup>/<sub>2</sub>",
        1.5: "<sup>+3</sup>/<sub>2</sub>",
        -1.5: "<sup>-3</sup>/<sub>2</sub>",
        2.5: "<sup>+5</sup>/<sub>2</sub>",
        -2.5: "<sup>-5</sup>/<sub>2</sub>",
    }
    for orbital, pdos, energy in projections.get_pdos():
        orbital_data = orbital.get_orbital_dict()
        kind_name = orbital_data["kind_name"]
        atom_position = [round(i, 2) for i in orbital_data["position"]]
        if atom_position not in list_positions:
            list_positions.append(atom_position)
        try:
            orbital_name = orbital.get_name_from_quantum_numbers(
                orbital_data["angular_momentum"], orbital_data["magnetic_number"]
            ).lower()
            if orbital_name in dict_html:
                orbital_name_plotly = dict_html[orbital_name]
            else:
                orbital_name_plotly = orbital_name

        except AttributeError:
            orbital_name = "j {j} l {l} m_j{m_j}".format(
                j=orbital_data["total_angular_momentum"],
                l=orbital_data["angular_momentum"],
                m_j=orbital_data["magnetic_number"],
            )
            orbital_name_plotly = "j={j} <i>l</i>={l} m<sub>j</sub>={m_j}".format(
                j=dict_html[orbital_data["total_angular_momentum"]],
                l=orbital_data["angular_momentum"],
                m_j=dict_html[orbital_data["magnetic_number"]],
            )
        if not selected_atoms:
            if group_tag == "atoms" and plot_tag == "total":
                key = r"{var}".format(var=atom_position)
            elif group_tag == "kinds" and plot_tag == "total":
                key = r"{var1}".format(var1=kind_name)
            elif group_tag == "atoms" and plot_tag == "orbital":
                key = r"{var1}<br>{var2}-{var3}".format(
                    var1=atom_position, var2=kind_name, var3=orbital_name_plotly
                )
            elif group_tag == "kinds" and plot_tag == "orbital":
                key = r"{var1}-{var2}".format(var1=kind_name, var2=orbital_name_plotly)
            else:
                key = None

            if key:
                if key in _pdos:
                    _pdos[key][1] += pdos
                else:
                    _pdos[key] = [energy, pdos]

        else:
            try:
                index = list_positions.index(atom_position)
                if index in selected_atoms:
                    if group_tag == "atoms" and plot_tag == "total":
                        key = r"{var}".format(var=atom_position)
                    elif group_tag == "kinds" and plot_tag == "total":
                        key = r"{var1}".format(var1=kind_name)
                    elif group_tag == "atoms" and plot_tag == "orbital":
                        key = r"{var1}<br>{var2}-{var3}".format(
                            var1=atom_position, var2=kind_name, var3=orbital_name_plotly
                        )
                    elif group_tag == "kinds" and plot_tag == "orbital":
                        key = r"{var1}-{var2}".format(
                            var1=kind_name, var2=orbital_name_plotly
                        )
                    else:
                        key = None

                    if key:
                        if key in _pdos:
                            _pdos[key][1] += pdos
                        else:
                            _pdos[key] = [energy, pdos]

            except ValueError:
                pass

    dos = []
    for label, (energy, pdos) in _pdos.items():
        if spin_type == "down":
            # invert y-axis
            pdos = -pdos
            label = f"{label} (↓)"

        if spin_type == "up":
            label = f"{label} (↑)"

        orbital_pdos = {
            "label": label,
            "x": energy.tolist(),
            "y": pdos.tolist(),
            "borderColor": cmap(label),
            "lineStyle": line_style,
        }
        dos.append(orbital_pdos)

    return dos
