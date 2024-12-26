"""Widgets for the QE app.

Authors: AiiDAlab team
"""

import base64
import hashlib
import os
import typing as t
from copy import deepcopy
from queue import Queue
from tempfile import NamedTemporaryFile
from threading import Event, Lock, Thread
from time import time

import ase
import ipywidgets as ipw
import numpy as np
import traitlets
from IPython.display import HTML, Javascript, clear_output, display

from aiida.orm import CalcJobNode, load_code, load_node
from aiida.orm import Data as orm_Data
from aiidalab_qe.common.mvc import Model
from aiidalab_widgets_base import (
    ComputationalResourcesWidget,
    StructureExamplesWidget,
    WizardAppWidgetStep,
)
from aiidalab_widgets_base.utils import (
    StatusHTML,
    list_to_string_range,
    string_range_to_list,
)

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

    def __init__(self, num_min_lines=10, max_output_height="200px", **kwargs):  # noqa: ARG002
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

        link_id = f"dl_{digest}"

        display(
            HTML(
                f"""
            <html>
            <body>
            <a id="{link_id}" download="{self.filename}" href="data:text/plain;base64,{payload}" download>
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
                return []

        elif "remote_folder" in calcjob.outputs:
            try:
                fn_out = calcjob.base.attributes.get("output_filename")
                self.filename = fn_out
                with NamedTemporaryFile() as tmpfile:
                    calcjob.outputs.remote_folder.getfile(fn_out, tmpfile.name)
                    return tmpfile.read().decode().splitlines()
            except OSError:
                return []
        else:
            return []

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

    def __init__(self, title="", **kwargs):
        self.title = title

        self._status_message = StatusHTML()
        self.atom_selection = ipw.Text(
            placeholder="e.g. 1..5 8 10",
            description="Index of atoms",
            value="",
            style={"description_width": "100px"},
            layout={"width": "initial"},
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
                ("3D (bulk systems)", "xyz"),
                ("2D (surfaces, slabs, ...)", "xy"),
                ("1D (wires)", "x"),
                ("0D (molecules)", "molecule"),
            ],
            value="xyz",
            layout={"width": "initial"},
        )
        self.select_periodicity = ipw.Button(
            description="Select",
            button_style="primary",
            layout={"width": "100px"},
        )
        self.scroll_note = ipw.HTML(
            value="<p style='font-style: italic;'>Note: The table is scrollable.</p>",
            layout={"visibility": "hidden"},
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
                    "<b>Set custom tags for atoms</b>",
                ),
                ipw.HTML(
                    """
                    <p>
                    These are used to distinguish atoms of the same chemical element. <br>
                    For example, they can be used to assign different initial magnetization values for antiferromagnetic systems.
                    </p>
                    <p style="font-weight: bold; color: #1f77b4;">NOTE:</p>
                    <ul style="padding-left: 2em; list-style-type: disc;">
                        <li>Atom indices start from 1, not 0. This means that the first atom in the list is numbered 1, the second atom is numbered 2, and so on.</li>
                    </ul>
                    </p>
                    """
                ),
                ipw.HBox(
                    [
                        self.atom_selection,
                        self.from_selection,
                        self.tag,
                    ]
                ),
                self.tag_display,
                self.scroll_note,
                ipw.HBox([self.add_tags, self.reset_tags, self.reset_all_tags]),
                self._status_message,
                ipw.HTML(
                    '<div style="margin-top: 20px;"><b>Set structure periodicity</b></div>'
                ),
                ipw.HTML("""
                    <p>Select the periodicity of your system.</p>
                    <p style="font-weight: bold; color: #1f77b4;">NOTE:</p>

                    <ul style="padding-left: 2em; list-style-type: disc;">
                        <li>For <b>2D</b> systems (e.g., surfaces, slabs), the non-periodic direction must be the third lattice vector (z-axis).</li>
                        <li>For <b>1D</b> systems (e.g., wires), the periodic direction must be the first lattice vector (x-axis).</li>

                    </ul>
                """),
                self.periodicity,
                self.select_periodicity,
            ],
            **kwargs,
        )

    def _display_table(self, _=None):
        """Function to control tag_display
        When given a list of atom in selection it will display a HTML table with Index, Element and Tag
        """
        selection = string_range_to_list(self.atom_selection.value)[0]
        current_tags = self.structure.get_tags()
        chemichal_symbols = self.structure.get_chemical_symbols()

        if (
            selection
            and (min(selection) >= 0)
            and (max(selection) <= (len(self.structure) - 1))
        ):
            table_data = []
            for index in selection:
                tag = current_tags[index]
                symbol = chemichal_symbols[index]
                if tag == 0:
                    tag = ""
                table_data.append([f"{index+ 1}", f"{symbol}", f"{tag}"])

            # Create an HTML table
            table_html = "<table>"
            table_html += "<tr><th>Index</th><th>Element</th><th>Tag</th></tr>"
            for row in table_data:
                table_html += "<tr>"
                for cell in row:
                    table_html += f"<td>{cell}</td>"
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
            self.scroll_note.layout = {"visibility": "visible"}
        else:
            self.tag_display.layout = {}
            with self.tag_display:
                clear_output()
            self.scroll_note.layout = {"visibility": "hidden"}

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
            "molecule": (False, False, False),
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

    def update_resources(self, change):
        if change["new"]:
            self.set_resource_defaults(load_code(change["new"]).computer)

    def set_resource_defaults(self, computer=None):
        if computer is None:
            self.num_nodes.disabled = True
            self.num_nodes.value = 1
            self.num_cpus.max = 1
            self.num_cpus.value = 1
            self.num_cpus.description = "CPUs"
        else:
            default_mpiprocs = computer.get_default_mpiprocs_per_machine()
            self.num_nodes.disabled = (
                True if computer.hostname == "localhost" else False
            )
            self.num_cpus.max = default_mpiprocs
            self.num_cpus.value = (
                1 if computer.hostname == "localhost" else default_mpiprocs
            )
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
        if "max_wallclock_seconds" in parameters:
            self.resource_detail.max_wallclock_seconds.value = parameters[
                "max_wallclock_seconds"
            ]

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
        Specify the parameters for the scheduler (only for advanced user). <br>
        These should be specified accordingly to the computer where the code will run.
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
        self.max_wallclock_seconds = ipw.BoundedIntText(
            value=3600 * 12,
            step=3600,
            min=60 * 10,
            max=3600 * 24,
            description="max seconds",
            style={"description_width": "100px"},
        )
        super().__init__(
            children=[
                self.prompt,
                self.ntasks_per_node,
                self.cpus_per_task,
                self.max_wallclock_seconds,
            ],
            **kwargs,
        )

    @property
    def parameters(self):
        return self.get_parameters()

    def get_parameters(self):
        """Return the parameters."""
        return {
            "ntasks_per_node": self.ntasks_per_node.value,
            "cpus_per_task": self.cpus_per_task.value,
            "max_wallclock_seconds": self.max_wallclock_seconds.value,
        }

    @parameters.setter
    def parameters(self, parameters):
        self.ntasks_per_node.value = parameters.get("ntasks_per_node", 1)
        self.cpus_per_task.value = parameters.get("cpus_per_task", 1)
        self.max_wallclock_seconds.value = parameters.get(
            "max_wallclock_seconds", 3600 * 12
        )

    def reset(self):
        """Reset the settings."""
        self.ntasks_per_node.value = 1
        self.cpus_per_task.value = 1
        self.max_wallclock_seconds.value = 3600 * 12


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
            ],
            **kwargs,
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


class LoadingWidget(ipw.HBox):
    """Widget for displaying a loading spinner."""

    def __init__(self, message="Loading", **kwargs):
        super().__init__(
            children=[
                ipw.Label(message),
                ipw.HTML("<i class='fa fa-spinner fa-spin fa-2x fa-fw'/>"),
            ],
            **kwargs,
        )
        self.add_class("loading")


class LazyLoader(ipw.VBox):
    identifier = "widget"

    def __init__(self, widget_class, widget_kwargs=None, **kwargs):
        super().__init__(
            children=[LoadingWidget(f"Loading {self.identifier}")],
            **kwargs,
        )

        self._widget_class = widget_class
        self._widget_kwargs = widget_kwargs or {}

        self.rendered = False

    def set_widget_kwargs(self, kwargs):
        self._widget_kwargs = kwargs

    def render(self):
        if self.rendered:
            return
        self.widget = self._widget_class(**self._widget_kwargs)
        self.children = [self.widget]
        self.rendered = True


class LazyLoadedStructureImporter(ipw.VBox):
    """A wrapper that may be used to lazy-load a structure import widget."""

    warning_message = "This may take some time to load"

    structure = traitlets.Union(
        [
            traitlets.Instance(ase.Atoms),
            traitlets.Instance(orm_Data),
        ],
        allow_none=True,
    )

    def __init__(self, title=None, **kwargs):
        self.title = title or "Structure importer"

        render_button = ipw.Button(
            layout=ipw.Layout(margin="0 10px 0 0", width="auto"),
            description=f"Load {self.title} widget",
            icon="refresh",
        )
        render_button.on_click(self._on_render_click)

        super().__init__(
            children=[
                ipw.HBox(
                    layout=ipw.Layout(align_items="center"),
                    children=[
                        render_button,
                        ipw.HTML("<span class='warning'>WARNING: </span>"),
                        ipw.HTML(f"<span>{self.warning_message}</span>"),
                    ],
                ),
            ],
            **kwargs,
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.children = [LoadingWidget(f"Loading {self.title} widget")]

        widget = self._get_widget()
        ipw.dlink(
            (widget, "structure"),
            (self, "structure"),
        )

        self.children = [widget]

        self.rendered = True

    def _get_widget(self):
        raise NotImplementedError()

    def _on_render_click(self, _):
        self.render()


class LazyLoadedOptimade(LazyLoadedStructureImporter):
    warning_message = (
        "OPTIMADE may take some time to load depending on your internet connection"
    )

    def _get_widget(self):
        from aiidalab_widgets_base.databases import OptimadeQueryWidget

        return OptimadeQueryWidget(embedded=False)


class LazyLoadedStructureBrowser(LazyLoadedStructureImporter):
    warning_message = (
        "The browser may take some time to load depending on the size of your database"
    )

    def _get_widget(self):
        from aiida import orm
        from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
        from aiidalab_widgets_base.structures import StructureBrowserWidget

        return StructureBrowserWidget(
            title=self.title,
            query_types=(
                orm.StructureData,
                orm.CifData,
                HubbardStructureData,
            ),
        )


class CategorizedStructureExamplesWidget(StructureExamplesWidget):
    """Extended widget to provide categorized example structures."""

    def __init__(self, examples_by_category, title="", **kwargs):
        self.examples_by_category = examples_by_category
        self._category_buttons = ipw.ToggleButtons(
            options=list(examples_by_category.keys()),
            description="Category:",
        )
        self._category_buttons.observe(self._on_category_change, names="value")

        # Initialize with the first category
        initial_category = next(iter(examples_by_category.keys()))
        super().__init__(
            examples=examples_by_category[initial_category], title=title, **kwargs
        )

        self.children = [self._category_buttons, *list(self.children)]

    def _on_category_change(self, change):
        """Update the dropdown when the category changes."""
        new_category = change["new"]
        new_examples = self.examples_by_category.get(new_category, [])
        self._select_structure.options = self.get_example_structures(new_examples)


class LinkButton(ipw.HTML):
    disabled = traitlets.Bool(False)

    def __init__(
        self,
        description=None,
        link="",
        in_place=False,
        classes="",
        icon="",
        disabled=False,
        **kwargs,
    ):
        super().__init__(**kwargs)

        html = f"""
            <a
                role="button"
                href="{link}"
                target="{"_self" if in_place else "_blank"}"
                class="{classes}"
            >
        """
        if icon:
            html += f"<i class='fa fa-{icon}'></i>"

        html += f"{description}</a>"

        self.value = html

        self.add_class("jupyter-button")
        self.add_class("link-button")

        self.disabled = disabled

    @traitlets.observe("disabled")
    def _on_disabled(self, change):
        if change["new"]:
            self.add_class("disabled")
        else:
            self.remove_class("disabled")


class QeWizardStepModel(Model):
    identifier = "QE wizard"


QWSM = t.TypeVar("QWSM", bound=QeWizardStepModel)


class QeWizardStep(ipw.VBox, WizardAppWidgetStep, t.Generic[QWSM]):
    def __init__(self, model: QWSM, **kwargs):
        self.loading_message = LoadingWidget(f"Loading {model.identifier} step")
        super().__init__(children=[self.loading_message], **kwargs)
        self._model = model
        self.rendered = False

    def render(self):
        if self.rendered:
            return
        self._render()
        self.rendered = True
        self._post_render()

    def _render(self):
        raise NotImplementedError()

    def _post_render(self):
        pass


class QeDependentWizardStep(QeWizardStep[QWSM]):
    missing_information_warning = "Missing information"

    previous_step_state = traitlets.UseEnum(WizardAppWidgetStep.State)

    def __init__(self, model: QWSM, **kwargs):
        super().__init__(model, **kwargs)
        self.previous_children = list(self.children)
        self.warning_message = ipw.HTML(
            f"""
            <div class="alert alert-danger">
                <b>Warning:</b> {self.missing_information_warning}
            </div>
        """
        )

    def render(self):
        if "PYTEST_CURRENT_TEST" in os.environ:
            super().render()
            return
        if self.previous_step_state is WizardAppWidgetStep.State.SUCCESS:
            self._hide_missing_information_warning()
            if not self.rendered:
                super().render()
                self.previous_children = list(self.children)
        else:
            self._show_missing_information_warning()

    def _show_missing_information_warning(self):
        self.children = [self.warning_message]
        self.rendered = False

    def _hide_missing_information_warning(self):
        self.children = self.previous_children
