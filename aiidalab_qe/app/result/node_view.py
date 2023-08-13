"""Results view widgets (MOVE TO OTHER MODULE!)

Authors: AiiDAlab team
"""

import shutil
import typing
from importlib import resources
from pathlib import Path
from tempfile import TemporaryDirectory

import ipywidgets as ipw
import traitlets
from aiida.cmdline.utils.common import get_workchain_report
from aiida.common import LinkType
from aiida.orm import CalcJobNode, WorkChainNode
from aiidalab_widgets_base import ProcessMonitor, register_viewer_widget
from filelock import FileLock, Timeout
from IPython.display import HTML, display
from jinja2 import Environment

from aiidalab_qe.app import static
from aiidalab_qe.app.plugins.relax.result import Result as RelaxResult
from aiidalab_qe.app.common.utils import get_entry_items

# I removed this `MinimalStructureViewer` class, because it is not used anywhere

# I removed all function related to pdos and bands to `pdos` and `bands` plugin folder.


class VBoxWithCaption(ipw.VBox):
    def __init__(self, caption, body, *args, **kwargs):
        super().__init__(children=[ipw.HTML(caption), body], *args, **kwargs)


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
        from .electronic_structure import ElectronicStructure
        from .summary_view import SummaryView

        if node.process_label != "QeAppWorkChain":
            super().__init__()
            return

        self.node = node
        ui_parameters = node.base.extras.get("ui_parameters", {})

        self.title = ipw.HTML(
            f"""
            <hr style="height:2px;background-color:#2097F3;">
            <h4>QE App Workflow (pk: {self.node.pk}) &mdash;
                {self.node.inputs.structure.get_formula()}
            </h4>
            """
        )
        self.workflows_summary = SummaryView(self.node)
        self.structure_tab = RelaxResult(self.node)
        self.electronic_tab = ElectronicStructure(self.node)

        self.summary_tab = ipw.VBox(children=[self.workflows_summary])
        self.result_tabs = ipw.Tab(
            children=[self.summary_tab, self.structure_tab, self.electronic_tab]
        )

        self.result_tabs.set_title(0, "Workflow Summary")
        self.result_tabs.set_title(1, "Final Geometry")
        self.result_tabs.set_title(2, "Electronic Structure")
        # add plugin specific settings
        entries = get_entry_items("aiidalab_qe.property", "result")
        self.results = {}
        for name, entry_point in entries.items():
            if not ui_parameters["workflow"]["properties"][name]:
                continue
            result = entry_point(self.node)
            self.results[name] = result
            self.result_tabs.children += (result,)
            self.result_tabs.set_title(len(self.result_tabs.children) - 1, result.title)

        # An ugly fix to the structure appearance problem
        # https://github.com/aiidalab/aiidalab-qe/issues/69
        def on_selected_index_change(change):
            index = change["new"]
            # Accessing the viewer only if the corresponding tab is present.
            if self.result_tabs._titles[str(index)] == "Final Geometry":
                if not hasattr(self.structure_tab, "_structure_view"):
                    return
                self.structure_tab._structure_view._viewer.handle_resize()

                def toggle_camera():
                    """Toggle camera between perspective and orthographic."""
                    self.structure_tab._structure_view._viewer.camera = (
                        "perspective"
                        if self.structure_tab._structure_view._viewer.camera
                        == "orthographic"
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
            self._show_workflow_output()
            # update the plugin specific results
            for result in self.result_tabs.children[1:]:
                # check if the result is already shown
                # check if the plugin workchain result is in the outputs
                if (
                    result.identifier not in self._results_shown
                    and result.identifier in self.node.outputs
                ):
                    result._update_view()
                    self._results_shown.add(result.identifier)
            # electronic structure needs to combine the results from the bands and pdos
            # TODO find a better way to show the results from several workchains
            self.result_tabs.children[2]._update_view()

    # I moved the _show_structure to the relax plugin part.
    # I moved the _show_electronic_structure to electronic_structure.py

    def _show_workflow_output(self):
        self.workflows_summary._update_view()
        if self.node.is_finished:
            self.workflows_output = WorkChainOutputs(self.node)
            self.result_tabs.children[0].children = [
                self.workflows_summary,
                self.workflows_output,
            ]
        else:
            self.result_tabs.children[0].children = [
                self.workflows_summary,
            ]
