"""Results view widgets (MOVE TO OTHER MODULE!)

Authors: AiiDAlab team
"""

import json
import random
import shutil
import typing
from importlib import resources
from pathlib import Path
from tempfile import TemporaryDirectory

import ipywidgets as ipw
import nglview
import traitlets
from aiida.cmdline.utils.common import get_workchain_report
from aiida.common import LinkType
from aiida.orm import CalcJobNode, Node, ProjectionData, WorkChainNode
from aiidalab_widgets_base import ProcessMonitor, register_viewer_widget
from aiidalab_widgets_base.viewers import StructureDataViewer
from ase import Atoms
from filelock import FileLock, Timeout
from IPython.display import HTML, display
from jinja2 import Environment
from monty.json import jsanitize
from traitlets import Instance, Int, List, Unicode, Union, default, observe, validate
from widget_bandsplot import BandsPlotWidget

from aiidalab_qe.app import static
from aiidalab_qe.app.report import generate_report_html


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
        self.wc_node = wc_node

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
            if (
                "structure" not in self._results_shown
                and "structure" in self.node.outputs
            ):
                self._show_structure()
                self._results_shown.add("structure")

            if "electronic_structure" not in self._results_shown and (
                "band_structure" in self.node.outputs or "dos" in self.node.outputs
            ):
                self._show_electronic_structure()
                self._results_shown.add("electronic_structure")

    def _show_structure(self):
        self._structure_view = StructureDataViewer(
            structure=self.node.outputs.structure
        )
        self.result_tabs.children[1].children = [self._structure_view]
        self.result_tabs.set_title(1, "Final Geometry")

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
