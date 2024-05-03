"""Widgets for the QE app.

Authors: AiiDAlab team
"""

import base64
import hashlib
from copy import deepcopy
from queue import Queue
from tempfile import NamedTemporaryFile
from threading import Event, Lock, Thread
from time import time

import ase
import ipywidgets as ipw
import numpy as np
import traitlets
from aiida.orm import CalcJobNode
from aiida.orm import Data as orm_Data
from aiida.orm import load_code, load_node
from aiidalab_widgets_base import ComputationalResourcesWidget
from aiidalab_widgets_base.utils import (
    StatusHTML,
    list_to_string_range,
    string_range_to_list,
)
from IPython.display import HTML, Javascript, clear_output, display
from pymatgen.core.periodic_table import Element
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


__all__ = [
    "CalcJobOutputFollower",
    "LogOutputWidget",
]


class RollingOutput(ipw.VBox):
    style = (
        "background-color: #253239; color: #cdd3df; line-height: normal; custom=test"
    )

    value = traitlets.Unicode()
    auto_scroll = traitlets.Bool()

    def __init__(self, num_min_lines=10, max_output_height="200px", **kwargs):
        self._num_min_lines = num_min_lines
        self._output = ipw.HTML(layout=ipw.Layout(min_width="50em"))
        self._refresh_output()
        super().__init__(
            [self._output],
            layout=ipw.Layout(max_height=max_output_height, min_width="51em"),
        )

    @traitlets.default("value")
    def _default_value(self):
        if self._num_min_lines > 0:
            return "\n" * self._num_min_lines

    @traitlets.default("auto_scroll")
    def _default_auto_scroll(self):
        return True

    def scroll_to_bottom(self):
        # Slight hack because it will scroll all widgets with the same class
        # name and max height to the bottom. That would primarily be an issue in
        # case that there are two LogOutputWidgets in the DOM. Could probably be
        # alleviated by adding a custom class.
        display(
            Javascript(
                """
            Array.from(document.getElementsByClassName('{class_name}'))
            .filter(el => el.style["max-height"] === "{max_height}")
            .forEach(el => el.scrollTop = el.scrollHeight)""".format(
                    class_name="p-Widget p-Panel jupyter-widgets widget-container widget-box widget-vbox",
                    max_height=self.layout.max_height,
                )
            )
        )

    @traitlets.observe("value")
    def _refresh_output(self, _=None):
        self._output.value = self._format_output(self.value)
        if self.auto_scroll:
            self.scroll_to_bottom()

    def _format_output(self, text):
        lines = text.splitlines()

        # Add empty lines to reach the minimum number of lines.
        lines += [""] * max(0, self._num_min_lines - len(lines))

        # Replace empty lines with single white space to ensure that they are
        # actually shown.
        lines = [line if len(line) > 0 else " " for line in lines]

        text = "\n".join(lines)
        return f"""<pre style="{self.style}">{text}</pre>"""


class DownloadButton(ipw.Button):
    # Adapted from https://stackoverflow.com/a/68683463
    """A Button widget for downloads with dynamic content."""

    filename = traitlets.Unicode()
    payload = traitlets.Bytes()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.on_click(self.__on_click)

    @traitlets.default("icon")
    def _default_icon(self):
        return "download"

    @traitlets.default("tooltip")
    def _default_tooltip(self):
        return "Download"

    def __on_click(self, _):
        digest = hashlib.md5(self.payload).hexdigest()  # bypass browser cache
        payload = base64.b64encode(self.payload).decode()

        id = f"dl_{digest}"

        display(
            HTML(
                f"""
            <html>
            <body>
            <a id="{id}" download="{self.filename}" href="data:text/plain;base64,{payload}" download>
            </a>

            <script>
            (function download() {{
            document.getElementById('{id}').click();
            }})()
            </script>

            </body>
            </html>
            """
            )
        )


class FilenameDisplayWidget(ipw.Box):
    value = traitlets.Unicode()

    def __init__(self, max_width=None, **kwargs):
        self.max_width = max_width
        self._html = ipw.HTML(
            layout={"margin": "0 0 0 50px"},
        )
        super().__init__([self._html], **kwargs)

    @traitlets.observe("value")
    def _observe_filename(self, change):
        icon = '<i class="fa fa-file-text-o" aria-hidden="true"></i>'
        width_style = f"width:{self.max_width};" if self.max_width else ""
        self._html.value = f"""
        <div style="
            white-space:nowrap;
            overflow:hidden;
            text-overflow:ellipsis;
            {width_style}">
            {icon} {change['new']}
        </div>
        """


class LogOutputWidget(ipw.VBox):
    filename = traitlets.Unicode()
    value = traitlets.Unicode()

    def __init__(self, placeholder=None, **kwargs):
        self.placeholder = placeholder

        self._rolling_output = RollingOutput(layout=ipw.Layout(flex="1 1 auto"))
        ipw.dlink(
            (self, "value"),
            (self._rolling_output, "value"),
            lambda value: value or self.placeholder or "",
        )

        self._filename_display = FilenameDisplayWidget(
            layout=ipw.Layout(width="auto"), max_width="55em"
        )
        ipw.dlink(
            (self, "filename"),
            (self._filename_display, "value"),
            lambda value: value or "[no filename]",
        )

        self._btn_download = DownloadButton(
            layout=ipw.Layout(width="30px", flex="5 1 auto"),
            disabled=True,
        )
        ipw.dlink((self, "filename"), (self._btn_download, "filename"))
        ipw.dlink(
            (self, "value"),
            (self._btn_download, "payload"),
            transform=lambda value: value.encode("utf-8"),
        )

        self._btn_scroll_down = ipw.Button(
            icon="angle-double-down",
            tooltip="Scroll to bottom",
            layout=ipw.Layout(width="30px", flex="1 1 auto"),
            disabled=True,
        )
        self._btn_scroll_down.on_click(
            lambda _: self._rolling_output.scroll_to_bottom()
        )

        self._btn_auto_scroll = ipw.ToggleButton(
            icon="magic",
            tooltip="Auto-scroll with content",
            value=True,
            layout=ipw.Layout(flex="1 1 auto", width="30px"),
        )
        ipw.link(
            (self._btn_auto_scroll, "value"), (self._rolling_output, "auto_scroll")
        )

        self._btns = ipw.VBox(
            [self._btn_download, self._btn_scroll_down, self._btn_auto_scroll],
            layout=ipw.Layout(min_width="36px", flex_flow="column"),
        )

        super().__init__(
            [
                self._filename_display,
                ipw.HBox([self._rolling_output, self._btns]),
            ],
            **kwargs,
        )

    @traitlets.default("placeholder")
    def _default_placeholder(self):
        return "[empty]"

    @traitlets.observe("value")
    def _observe_value(self, change):
        self._btn_download.disabled = not change["new"]
        self._btn_scroll_down.disabled = not change["new"]


class CalcJobOutputFollower(traitlets.HasTraits):
    calcjob_uuid = traitlets.Unicode(allow_none=True)
    filename = traitlets.Unicode(allow_none=True)
    output = traitlets.List(trait=traitlets.Unicode)
    lineno = traitlets.Int()

    def __init__(self, **kwargs):
        self._output_queue = Queue()

        self._lock = Lock()
        self._push_thread = None
        self._pull_thread = None
        self._stop_follow_output = Event()
        self._follow_output_thread = None

        super().__init__(**kwargs)

    @traitlets.observe("calcjob_uuid")
    def _observe_calcjob(self, change):
        calcjob_uuid = change["new"]
        if change["old"] == calcjob_uuid:
            return

        with self._lock:
            # Stop following
            self._stop_follow_output.set()

            if self._follow_output_thread:
                self._follow_output_thread.join()
                self._follow_output_thread = None

            # Reset all traitlets and signals.
            self.output.clear()
            self.lineno = 0
            self._stop_follow_output.clear()

            # (Re/)start following
            if change["new"]:
                self._follow_output_thread = Thread(
                    target=self._follow_output, args=(calcjob_uuid,)
                )
                self._follow_output_thread.start()

    def _follow_output(self, calcjob_uuid):
        """Monitor calcjob and orchestrate pushing and pulling of output."""
        self._pull_thread = Thread(target=self._pull_output)
        self._pull_thread.start()
        self._push_thread = Thread(target=self._push_output, args=(calcjob_uuid,))
        self._push_thread.start()

    def _fetch_output(self, calcjob):
        assert isinstance(calcjob, CalcJobNode)
        if "retrieved" in calcjob.outputs:
            try:
                self.filename = calcjob.base.attributes.get("output_filename")
                with calcjob.outputs.retrieved.base.repository.open(self.filename) as f:
                    return f.read().splitlines()
            except OSError:
                return list()

        elif "remote_folder" in calcjob.outputs:
            try:
                fn_out = calcjob.base.attributes.get("output_filename")
                self.filename = fn_out
                with NamedTemporaryFile() as tmpfile:
                    calcjob.outputs.remote_folder.getfile(fn_out, tmpfile.name)
                    return tmpfile.read().decode().splitlines()
            except OSError:
                return list()
        else:
            return list()

    _EOF = None

    def _push_output(self, calcjob_uuid, delay=0.2):
        """Push new log lines onto the queue."""
        lineno = 0
        calcjob = load_node(calcjob_uuid)
        while True:
            try:
                lines = self._fetch_output(calcjob)
            except Exception as error:
                self._output_queue.put([f"[ERROR: {error}]"])
            else:
                self._output_queue.put(lines[lineno:])
                lineno = len(lines)
            finally:
                if calcjob.is_sealed or self._stop_follow_output.wait(delay):
                    # Pushing EOF signals to the pull thread to stop.
                    self._output_queue.put(self._EOF)
                    break  # noqa: B012

    def _pull_output(self):
        """Pull new log lines from the queue and update traitlets."""
        while True:
            item = self._output_queue.get()
            if item is self._EOF:
                self._output_queue.task_done()
                break
            else:  # item is 'new lines'
                with self.hold_trait_notifications():
                    self.output.extend(item)
                    self.lineno += len(item)
                self._output_queue.task_done()


class ProgressBar(ipw.HBox):
    class AnimationRate(float):
        pass

    description = traitlets.Unicode()
    value = traitlets.Union([traitlets.Float(), traitlets.Instance(AnimationRate)])
    bar_style = traitlets.Unicode()

    _animation_rate = traitlets.Float()

    def __init__(self, description_layout=None, *args, **kwargs):
        if description_layout is None:
            description_layout = ipw.Layout(width="auto", flex="2 1 auto")

        self._label = ipw.Label(layout=description_layout)
        self._progress_bar = ipw.FloatProgress(
            min=0, max=1.0, layout=ipw.Layout(width="auto", flex="1 1 auto")
        )

        traitlets.link((self, "description"), (self._label, "value"))
        traitlets.link((self, "bar_style"), (self._progress_bar, "bar_style"))

        self._animate_stop_event = Event()
        self._animate_thread = None

        super().__init__([self._label, self._progress_bar], *args, **kwargs)

    def _animate(self, refresh_rate=0.01):
        v0 = self._progress_bar.value
        t0 = time()

        while not self._animate_stop_event.wait(refresh_rate):
            self._progress_bar.value = (v0 + (time() - t0) * self._animation_rate) % 1.0

    def _start_animate(self):
        if self._animate_thread is not None:
            raise RuntimeError("Cannot start animation more than once!")

        self._animate_thread = Thread(target=self._animate)
        self._animate_thread.start()

    def _stop_animate(self):
        self._animate_stop_event.set()
        self._animate_thread.join()
        self._animate_stop_event.clear()
        self._animate_thread = None

    @traitlets.default("_animation_rate")
    def _default_animation_rate(self):
        return 0

    @traitlets.observe("_animation_rate")
    def _observe_animation_rate(self, change):
        if change["new"] and not change["old"]:
            self._start_animate()
        elif not change["new"] and change["old"]:
            self._stop_animate()

    @traitlets.validate("value")
    def _validate_value(self, proposal):
        if isinstance(proposal["value"], self.AnimationRate):
            if proposal["value"] < 0:
                raise traitlets.TraitError("The animation rate must be non-negative.")

        elif not 0 <= proposal["value"] <= 1.0:
            raise traitlets.TraitError("The value must be between 0 and 1.0.")

        return proposal["value"]

    @traitlets.observe("value")
    def _observe_value(self, change):
        if isinstance(change["new"], self.AnimationRate):
            self._animation_rate = change["new"]
        else:
            self._animation_rate = 0
            self._progress_bar.value = change["new"]


class AddingTagsEditor(ipw.VBox):
    """Editor for adding tags to atoms."""

    structure = traitlets.Instance(ase.Atoms, allow_none=True)
    selection = traitlets.List(traitlets.Int, allow_none=True)
    input_selection = traitlets.List(traitlets.Int, allow_none=True)
    structure_node = traitlets.Instance(orm_Data, allow_none=True, read_only=True)

    def __init__(self, title=""):
        self.title = title
        self._status_message = StatusHTML()
        self.atom_selection = ipw.Text(
            description="Index of atoms", value="", layout={"width": "initial"}
        )
        self.from_selection = ipw.Button(description="From selection")
        self.from_selection.on_click(self._from_selection)
        self.tag = ipw.BoundedIntText(
            description="Tag", value=1, min=0, max=4, layout={"width": "initial"}
        )
        self.add_tags = ipw.Button(
            description="Update tags",
            button_style="primary",
            layout={"width": "initial"},
        )

        self.reset_tags = ipw.Button(
            description="Reset tags",
            button_style="primary",
            layout={"width": "initial"},
        )
        self.reset_all_tags = ipw.Button(
            description="Reset all tags",
            button_style="warning",
            layout={"width": "initial"},
        )
        self.periodicity = ipw.RadioButtons(
            options=[
                "xyz",
                "xy",
                "x",
            ],
            value="xyz",
            description="Periodicty: ",
            layout={"width": "initial"},
        )
        self.select_periodicity = ipw.Button(
            description="Select",
            button_style="primary",
            layout={"width": "100px"},
        )
        self.tag_display = ipw.Output()
        self.add_tags.on_click(self._add_tags)
        self.reset_tags.on_click(self._reset_tags)
        self.reset_all_tags.on_click(self._reset_all_tags)
        self.atom_selection.observe(self._display_table, "value")
        self.add_tags.on_click(self._display_table)
        self.reset_tags.on_click(self._display_table)
        self.reset_all_tags.on_click(self._display_table)
        self.select_periodicity.on_click(self._select_periodicity)
        super().__init__(
            children=[
                ipw.HTML(
                    "<b>Adding a tag to atoms</b>",
                ),
                ipw.HBox(
                    [
                        self.atom_selection,
                        self.from_selection,
                        self.tag,
                    ]
                ),
                self.tag_display,
                ipw.HBox([self.add_tags, self.reset_tags, self.reset_all_tags]),
                self._status_message,
                ipw.HTML(
                    "<b>Define periodicity</b>",
                ),
                self.periodicity,
                self.select_periodicity,
            ]
        )

    def _display_table(self, _=None):
        """Function to control tag_display
        When given a list of atom in selection it will display a HTML table with Index, Element and Tag
        """
        selection = string_range_to_list(self.atom_selection.value)[0]
        current_tags = self.structure.get_tags()
        chemichal_symbols = self.structure.get_chemical_symbols()

        if selection and (max(selection) <= (len(self.structure) - 1)):
            table_data = []
            for index in selection:
                tag = current_tags[index]
                symbol = chemichal_symbols[index]
                if tag == 0:
                    tag = ""
                table_data.append(
                    ["{}".format(index), "{}".format(symbol), "{}".format(tag)]
                )

            # Create an HTML table
            table_html = "<table>"
            table_html += "<tr><th>Index</th><th>Element</th><th>Tag</th></tr>"
            for row in table_data:
                table_html += "<tr>"
                for cell in row:
                    table_html += "<td>{}</td>".format(cell)
                table_html += "</tr>"
            table_html += "</table>"

            # Set layout to a fix size
            self.tag_display.layout = {
                "overflow": "auto",
                "height": "100px",
                "width": "150px",
            }
            with self.tag_display:
                clear_output()
                display(HTML(table_html))
        else:
            self.tag_display.layout = {}
            with self.tag_display:
                clear_output()

    def _from_selection(self, _=None):
        """Set the atom selection from the current selection."""
        self.atom_selection.value = list_to_string_range(self.selection)

    def _add_tags(self, _=None):
        """Add tags to the selected atoms."""
        if not self.atom_selection.value:
            self._status_message.message = """
            <div class="alert alert-info">
            <strong>Please select atoms first.</strong>
            </div>
            """
        else:
            selection = string_range_to_list(self.atom_selection.value)[0]
            new_structure = deepcopy(self.structure)
            if not new_structure.get_tags().any():
                new_tags = np.zeros(len(new_structure))
            else:
                new_tags = new_structure.get_tags()
            new_tags[selection] = self.tag.value
            new_structure.set_tags(new_tags)
            self.structure = None
            self.structure = deepcopy(new_structure)
            self.input_selection = None
            self.input_selection = deepcopy(self.selection)

    def _reset_tags(self, _=None):
        """Clear tags from selected atoms."""
        if not self.atom_selection.value:
            self._status_message.message = """
            <div class="alert alert-info">
            <strong>Please select atoms first.</strong>
            </div>
            """
        else:
            selection = string_range_to_list(self.atom_selection.value)[0]
            new_structure = deepcopy(self.structure)
            new_tags = new_structure.get_tags()
            new_tags[selection] = 0
            new_structure.set_tags(new_tags)
            self.structure = None
            self.structure = deepcopy(new_structure)
            self.input_selection = None
            self.input_selection = deepcopy(self.selection)

    def _reset_all_tags(self, _=None):
        """Clear all tags."""
        new_structure = deepcopy(self.structure)
        new_tags = np.zeros(len(new_structure))
        new_structure.set_tags(new_tags)
        self.structure = None
        self.structure = deepcopy(new_structure)
        self.input_selection = None
        self.input_selection = deepcopy(self.selection)

    def _select_periodicity(self, _=None):
        """Select periodicity."""
        periodicity_options = {
            "xyz": (True, True, True),
            "xy": (True, True, False),
            "x": (True, False, False),
        }
        new_structure = deepcopy(self.structure)
        new_structure.set_pbc(periodicity_options[self.periodicity.value])
        self.structure = None
        self.structure = deepcopy(new_structure)


class QEAppComputationalResourcesWidget(ipw.VBox):
    value = traitlets.Unicode(allow_none=True)
    nodes = traitlets.Int(default_value=1)
    cpus = traitlets.Int(default_value=1)

    def __init__(self, **kwargs):
        """Widget to setup the compute resources, which include the code,
        the number of nodes and the number of cpus.
        """
        self.code_selection = ComputationalResourcesWidget(**kwargs)
        self.code_selection.layout.width = "80%"

        self.num_nodes = ipw.BoundedIntText(
            value=1, step=1, min=1, max=1000, description="Nodes", width="10%"
        )
        self.num_cpus = ipw.BoundedIntText(
            value=1, step=1, min=1, description="CPUs", width="10%"
        )
        self.btn_setup_resource_detail = ipw.ToggleButton(description="More")
        self.btn_setup_resource_detail.observe(self._setup_resource_detail, "value")
        self._setup_resource_detail_output = ipw.Output(layout={"width": "500px"})

        # combine code, nodes and cpus
        children = [
            ipw.HBox(
                children=[
                    self.code_selection,
                    self.num_nodes,
                    self.num_cpus,
                    self.btn_setup_resource_detail,
                ]
            ),
            self._setup_resource_detail_output,
        ]
        super().__init__(children=children, **kwargs)

        self.resource_detail = ResourceDetailSettings()
        traitlets.dlink(
            (self.num_cpus, "value"), (self.resource_detail.ntasks_per_node, "value")
        )
        traitlets.link((self.code_selection, "value"), (self, "value"))

    @traitlets.observe("value")
    def _update_resources(self, change):
        if change["new"]:
            self.set_resource_defaults(load_code(change["new"]).computer)

    def set_resource_defaults(self, computer=None):
        import os

        if computer is None or computer.hostname == "localhost":
            self.num_nodes.disabled = True
            self.num_nodes.value = 1
            self.num_cpus.max = os.cpu_count()
            self.num_cpus.value = 1
            self.num_cpus.description = "CPUs"
        else:
            default_mpiprocs = computer.get_default_mpiprocs_per_machine()
            self.num_nodes.disabled = False
            self.num_cpus.max = default_mpiprocs
            self.num_cpus.value = default_mpiprocs
            self.num_cpus.description = "CPUs"

    @property
    def parameters(self):
        return self.get_parameters()

    def get_parameters(self):
        """Return the parameters."""
        parameters = {
            "code": self.code_selection.value,
            "nodes": self.num_nodes.value,
            "cpus": self.num_cpus.value,
        }
        parameters.update(self.resource_detail.parameters)
        return parameters

    @parameters.setter
    def parameters(self, parameters):
        self.set_parameters(parameters)

    def set_parameters(self, parameters):
        """Set the parameters."""
        self.code_selection.value = parameters["code"]
        if "nodes" in parameters:
            self.num_nodes.value = parameters["nodes"]
        if "cpus" in parameters:
            self.num_cpus.value = parameters["cpus"]
        if "ntasks_per_node" in parameters:
            self.resource_detail.ntasks_per_node.value = parameters["ntasks_per_node"]
        if "cpus_per_task" in parameters:
            self.resource_detail.cpus_per_task.value = parameters["cpus_per_task"]

    def _setup_resource_detail(self, _=None):
        with self._setup_resource_detail_output:
            clear_output()
            if self.btn_setup_resource_detail.value:
                self._setup_resource_detail_output.layout = {
                    "width": "500px",
                    "border": "1px solid gray",
                }

                children = [
                    self.resource_detail,
                ]
                display(*children)
            else:
                self._setup_resource_detail_output.layout = {
                    "width": "500px",
                    "border": "none",
                }


class ResourceDetailSettings(ipw.VBox):
    """Widget for setting the Resource detail."""

    prompt = ipw.HTML(
        """<div style="line-height:120%; padding-top:0px">
        <p style="padding-bottom:10px">
        Specify the parameters for the scheduler (only for advanced user).
        </p></div>"""
    )

    def __init__(self, **kwargs):
        self.ntasks_per_node = ipw.BoundedIntText(
            value=1,
            step=1,
            min=1,
            max=1000,
            description="ntasks-per-node",
            style={"description_width": "100px"},
        )
        self.cpus_per_task = ipw.BoundedIntText(
            value=1,
            step=1,
            min=1,
            description="cpus-per-task",
            style={"description_width": "100px"},
        )
        super().__init__(
            children=[self.prompt, self.ntasks_per_node, self.cpus_per_task], **kwargs
        )

    @property
    def parameters(self):
        return self.get_parameters()

    def get_parameters(self):
        """Return the parameters."""
        return {
            "ntasks_per_node": self.ntasks_per_node.value,
            "cpus_per_task": self.cpus_per_task.value,
        }

    @parameters.setter
    def parameters(self, parameters):
        self.ntasks_per_node.value = parameters.get("ntasks_per_node", 1)
        self.cpus_per_task.value = parameters.get("cpus_per_task", 1)

    def reset(self):
        """Reset the settings."""
        self.ntasks_per_node.value = 1
        self.cpus_per_task.value = 1


class ParallelizationSettings(ipw.VBox):
    """Widget for setting the parallelization settings."""

    prompt = ipw.HTML(
        """<div style="line-height:120%; padding-top:0px">
        <p style="padding-bottom:10px">
        Specify the number of k-points pools for the pw.x calculations (only for advanced user).
        </p></div>"""
    )

    def __init__(self, **kwargs):
        extra = {
            "style": {"description_width": "150px"},
            "layout": {"min_width": "180px"},
        }
        self.npool = ipw.BoundedIntText(
            value=1, step=1, min=1, max=128, description="Number of k-pools", **extra
        )
        self.override = ipw.Checkbox(
            escription="",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="20px"),
        )
        self.override.observe(self.set_visibility, "value")
        super().__init__(
            children=[
                ipw.HBox(
                    children=[self.override, self.prompt, self.npool],
                    layout=ipw.Layout(justify_content="flex-start"),
                ),
            ]
        )
        # set the default visibility of the widget
        self.npool.layout.display = "none"

    def set_visibility(self, change):
        if change["new"]:
            self.npool.layout.display = "block"
        else:
            self.npool.layout.display = "none"

    def reset(self):
        """Reset the parallelization settings."""
        self.npool.value = 1


class PwCodeResourceSetupWidget(QEAppComputationalResourcesWidget):
    """ComputationalResources Widget for the pw.x calculation."""

    nodes = traitlets.Int(default_value=1)

    def __init__(self, **kwargs):
        # By definition, npool must be a divisor of the total number of k-points
        # thus we can not set a default value here, or from the computer.
        self.parallelization = ParallelizationSettings()
        super().__init__(**kwargs)
        # add nodes and cpus into the children of the widget
        self.children += (self.parallelization,)

    def get_parallelization(self):
        """Return the parallelization settings."""
        parallelization = (
            {"npool": self.parallelization.npool.value}
            if self.parallelization.override.value
            else {}
        )
        return parallelization

    def set_parallelization(self, parallelization):
        """Set the parallelization settings."""
        if "npool" in parallelization:
            self.parallelization.override.value = True
            self.parallelization.npool.value = parallelization["npool"]

    def get_parameters(self):
        """Return the parameters."""
        parameters = super().get_parameters()
        parameters.update({"parallelization": self.get_parallelization()})
        return parameters

    def set_parameters(self, parameters):
        """Set the parameters."""
        super().set_parameters(parameters)
        if "parallelization" in parameters:
            self.set_parallelization(parameters["parallelization"])


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
                [ipw.HTML("Define U value [eV] ")]
                + widgets_list
                + [self.eigenvalues_help, self.eigenvalues_label]
            )
        else:
            hubbard_widget = ipw.VBox([ipw.HTML("Define U value [eV] ")] + widgets_list)
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
                        )
                    ]
                    + widgets_list_up,
                )

                row_down = ipw.HBox(
                    children=[
                        ipw.Label(
                            "Down:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="50px"
                            ),
                        )
                    ]
                    + widgets_list_down,
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
            parameters = self.set_parameters_from_hubbardstructe(self.input_structure)
            self.set_hubbard_widget(parameters)
            self.activate_hubbard.value = True

    def set_parameters_from_hubbardstructe(self, hubbard_structure):
        hubbard_parameters = hubbard_structure.hubbard.dict()["parameters"]
        parameters = {f"{hubbard_structure.sites[item['atom_index']].kind_name} - {item['atom_manifold']}": item['value'] for item in hubbard_parameters}
        return parameters

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
        self.input_structure = None
        self.hubbard_widget = self.create_hubbard_widget()
        self.eigen_values_widget = self.create_eigenvalues_widget()
        with self.hubbard_widget_out:
            clear_output()
        with self.eigen_values_widget_out:
            clear_output()

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
