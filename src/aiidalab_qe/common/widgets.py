"""Widgets for the QE app.

Authors: AiiDAlab team
"""

import base64
import hashlib
import subprocess
import typing as t
from copy import deepcopy
from queue import Queue
from tempfile import NamedTemporaryFile
from threading import Event, Lock, Thread
from time import time

import anywidget
import ase
import ipywidgets as ipw
import numpy as np
import requests as req
import traitlets
from IPython.display import HTML, Javascript, clear_output, display
from pymatgen.io.ase import AseAtomsAdaptor
from shakenbreak.distortions import distort, local_mc_rattle, rattle

from aiida.orm import CalcJobNode, load_code, load_node
from aiida.orm import Data as orm_Data
from aiidalab_widgets_base import (
    ComputationalResourcesWidget,
    StructureExamplesWidget,
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
        self._output = ipw.HTML(layout=ipw.Layout(min_width="80em"))
        self._refresh_output()
        super().__init__(
            children=[
                self._output,
            ],
            layout=ipw.Layout(max_height=max_output_height),
        )
        self.add_class("rolling-output")

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
        self._html = ipw.HTML()
        super().__init__([self._html], **kwargs)
        self.add_class("filename-display")

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
            {icon} {change["new"]}
        </div>
        """


class LogOutputWidget(ipw.VBox):
    filename = traitlets.Unicode()
    value = traitlets.Unicode()

    def __init__(self, placeholder=None, **kwargs):
        self.placeholder = placeholder

        self._rolling_output = RollingOutput(
            layout=ipw.Layout(flex="1 1 auto"),
            max_output_height="unset",
        )
        ipw.dlink(
            (self, "value"),
            (self._rolling_output, "value"),
            lambda value: value or self.placeholder or "",
        )

        self._filename_display = FilenameDisplayWidget(layout=ipw.Layout(width="auto"))
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
            children=[
                self._filename_display,
                ipw.HBox(
                    children=[
                        self._rolling_output,
                        self._btns,
                    ],
                    layout=ipw.Layout(height="100%"),
                ),
            ],
            layout=ipw.Layout(height="100%"),
            **kwargs,
        )
        self.add_class("log-output")

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
            description="Tag", value=1, min=0, max=11, layout={"width": "initial"}
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

        super().__init__(
            children=[
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
                table_data.append([f"{index + 1}", f"{symbol}", f"{tag}"])

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


class PeriodicityEditor(ipw.VBox):
    """Editor for changing periodicity of structures."""

    structure = traitlets.Instance(ase.Atoms, allow_none=True)

    def __init__(self, title="", **kwargs):
        self.title = title

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
        self.apply_periodicity = ipw.Button(
            description="Apply",
            button_style="primary",
            layout={"width": "100px"},
        )
        self.apply_periodicity.on_click(self._select_periodicity)

        super().__init__(
            children=[
                ipw.HTML("""
                    <p>Select the periodicity of your system.</p>
                    <p style="font-weight: bold; color: #1f77b4;">NOTE:</p>

                    <ul style="padding-left: 2em; list-style-type: disc;">
                        <li>For <b>2D</b> systems (e.g., surfaces, slabs), the non-periodic direction must be the third lattice vector (z-axis).</li>
                        <li>For <b>1D</b> systems (e.g., wires), the periodic direction must be the first lattice vector (x-axis).</li>

                    </ul>
                """),
                self.periodicity,
                self.apply_periodicity,
            ],
            **kwargs,
        )

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
        self.code_selection = ComputationalResourcesWidget(
            include_setup_widget=False,
            fetch_codes=True,  # TODO resolve testing issues when set to `False`
            **kwargs,
        )
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
        self.message = ipw.Label(message)
        super().__init__(
            children=[
                ipw.Label(message),
                ipw.HTML(
                    value="<i class='fa fa-spinner fa-spin fa-2x fa-fw'/>",
                    layout=ipw.Layout(margin="12px 0 6px"),
                ),
            ],
            layout=ipw.Layout(
                justify_content="center",
                align_items="center",
                **kwargs.pop("layout", {}),
            ),
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
        class_="",
        style_="",
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
                style="cursor: default; {style_}"
            >
        """
        if icon:
            html += f"<i class='fa fa-{icon}'></i>"

        html += f"{description}</a>"

        self.value = html

        self.add_class("jupyter-button")
        self.add_class("widget-button")
        self.add_class("link-button")
        self.add_class(class_)

        self.disabled = disabled

    @traitlets.observe("disabled")
    def _on_disabled(self, change):
        if change["new"]:
            self.add_class("disabled")
        else:
            self.remove_class("disabled")


class TableWidget(anywidget.AnyWidget):
    _esm = """
    function render({ model, el }) {
        let domElement = document.createElement("div");
        el.classList.add("custom-table");
        let selectedIndices = [];

        function drawTable() {
            const data = model.get("data");
            domElement.innerHTML = "";
            let innerHTML = '<table><tr>' + data[0].map(header => `<th>${header}</th>`).join('') + '</tr>';

            for (let i = 1; i < data.length; i++) {
                innerHTML += '<tr>' + data[i].map(cell => `<td>${cell}</td>`).join('') + '</tr>';
            }

            innerHTML += "</table>";
            domElement.innerHTML = innerHTML;

            const rows = domElement.querySelectorAll('tr');
            rows.forEach((row, index) => {
                if (index > 0) {
                    row.addEventListener('click', () => {
                        const rowIndex = index - 1;
                        if (selectedIndices.includes(rowIndex)) {
                            selectedIndices = selectedIndices.filter(i => i !== rowIndex);
                            row.classList.remove('selected-row');
                        } else {
                            selectedIndices.push(rowIndex);
                            row.classList.add('selected-row');
                        }
                        model.set('selected_rows', [...selectedIndices]);
                        model.save_changes();
                    });

                    row.addEventListener('mouseover', () => {
                        if (!row.classList.contains('selected-row')) {
                            row.classList.add('hover-row');
                        }
                    });

                    row.addEventListener('mouseout', () => {
                        row.classList.remove('hover-row');
                    });
                }
            });
        }

        function updateSelection() {
            const newSelection = model.get("selected_rows");
            selectedIndices = [...newSelection]; // Synchronize the JavaScript state with the Python state
            const rows = domElement.querySelectorAll('tr');
            rows.forEach((row, index) => {
                if (index > 0) {
                    if (selectedIndices.includes(index - 1)) {
                        row.classList.add('selected-row');
                    } else {
                        row.classList.remove('selected-row');
                    }
                }
            });
        }

        drawTable();
        model.on("change:data", drawTable);
        model.on("change:selected_rows", updateSelection);
        el.appendChild(domElement);
    }
    export default { render };
    """
    _css = """
    .custom-table table, .custom-table th, .custom-table td {
        border: 1px solid black;
        border-collapse: collapse;
        text-align: left;
        padding: 4px;
    }
    .custom-table th, .custom-table td {
        min-width: 50px;
        word-wrap: break-word;
    }
    .custom-table table {
        width: 70%;
        font-size: 1.0em;
    }
    /* Hover effect with light gray background */
    .custom-table tr.hover-row:not(.selected-row) {
        background-color: #f0f0f0;
    }
    /* Selected row effect with light green background */
    .custom-table tr.selected-row {
        background-color: #DFF0D8;
    }
    """
    data = traitlets.List().tag(sync=True)
    selected_rows = traitlets.List().tag(sync=True)


class HBoxWithUnits(ipw.HBox):
    def __init__(self, widget: ipw.ValueWidget, units: str, **kwargs):
        super().__init__(
            children=[
                widget,
                ipw.HTML(units),
            ],
            layout=ipw.Layout(
                align_items="center",
                grid_gap="2px",
            ),
            **kwargs,
        )


class ArchiveImporter(ipw.VBox):
    GITHUB = "https://github.com"
    INFO_TEMPLATE = "{} <i class='fa fa-spinner fa-spin'></i>"

    def __init__(
        self,
        repo: str,
        tag: str,
        archive_list: str,
        archives_dir: str,
        logger: t.Optional[dict] = None,
        **kwargs,
    ):
        refs = "refs/tags"
        raw = f"https://raw.githubusercontent.com/{repo}/{refs}/{tag}"
        self.archive_list_url = f"{raw}/{archive_list}"
        self.archives_url = f"{self.GITHUB}/{repo}/raw/{refs}/{tag}/{archives_dir}"
        if logger:
            self.logger_placeholder = logger.get("placeholder", "")
            self.clear_log_on_import = logger.get("clear_on_import", False)
            self.logger = RollingOutput()
            self.logger.value = self.logger_placeholder
        super().__init__(children=[LoadingWidget()], **kwargs)

    def render(self):
        self.selector = ipw.SelectMultiple(
            options=[],
            description="Examples:",
            rows=10,
            style={"description_width": "initial"},
            layout=ipw.Layout(width="auto"),
        )

        self.import_button = ipw.Button(
            description="Import",
            button_style="success",
            layout=ipw.Layout(width="fit-content"),
            icon="download",
        )
        ipw.dlink(
            (self.selector, "value"),
            (self.import_button, "disabled"),
            lambda value: not value,
        )
        self.import_button.on_click(self.import_archives)

        self.info = ipw.HTML()

        history_link = LinkButton(
            description="Calculation history",
            link="./calculation_history.ipynb",
            icon="list",
            class_="mod-primary",
            style_="color: white;",
            layout=ipw.Layout(
                width="fit-content",
                margin="2px 0 2px auto",
            ),
        )

        accordion = None
        if self.logger:
            accordion = ipw.Accordion(children=[self.logger])
            accordion.set_title(0, "Archive import log")
            accordion.selected_index = None

        self.children = [
            self.selector,
            ipw.HBox(
                [
                    self.import_button,
                    self.info,
                    history_link,
                ],
                layout=ipw.Layout(
                    margin="2px 2px 4px 68px",
                    align_items="center",
                    grid_gap="4px",
                ),
            ),
            accordion or ipw.Box,
        ]

        self.selector.options = self.get_options()

    def get_options(self) -> list[tuple[str, str]]:
        response: req.Response = req.get(self.archive_list_url)
        if not response.ok:
            self.info.value = "Failed to fetch archive list"
            return []
        if archives := response.content.decode("utf-8").strip().split("\n"):
            return [(archive, archive.split("-")[0].strip()) for archive in archives]
        self.info.value = "No archives found"
        return []

    def import_archives(self, _):
        self.import_button.disabled = True
        if self.logger and self.clear_log_on_import:
            self.logger.value = ""
        for filename in self.selector.value:
            self.import_archive(filename)
        self.import_button.disabled = False

    def import_archive(self, filename: str):
        self.info.value = self.INFO_TEMPLATE.format(f"Importing {filename}")
        file_url = f"{self.archives_url}/{filename}"
        process = subprocess.Popen(
            ["verdi", "archive", "import", "-v", "critical", file_url],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        stdout, stderr = process.communicate()
        self._report(filename, stdout, stderr)

    def _report(self, filename: str, stdout: str, stderr: str):
        if stderr and "Success" not in stdout:
            self.info.value = f"Failed to import {filename}"
        if self.logger:
            self.logger.value += stdout
            if stderr:
                if "Success" not in stdout:
                    self.logger.value += "\n[ERROR]"
                    self.info.value = "Error -> see log for details"
                else:
                    self.info.value = ""
                self.logger.value += f"\n{stderr}"


class ShakeNBreakEditor(ipw.VBox):
    structure = traitlets.Instance(ase.Atoms, allow_none=True)
    selection = traitlets.List(traitlets.Int)
    structure_node = traitlets.Instance(orm_Data, allow_none=True, read_only=True)

    def __init__(self, title="Editor ShakeNbreak"):
        self.title = title
        self._status_message = StatusHTML()
        self.num_nearest_neighbours = ipw.IntText(
            description="Num. neighbours atoms:",
            value=8,
            step=1,
            style={"description_width": "150px"},
        )
        self.defect_position = ipw.Text(
            description="Defect atom:",
            style={"description_width": "150px"},
        )
        # Only select one atom!
        self.btn_defect_position = ipw.Button(
            description="From selection",
            style={"description_width": "150px"},
        )
        # Center of the selection
        self.btn_defect_position_vac = ipw.Button(
            description="From selection",
            style={"description_width": "150px"},
        )
        self.vacancy_coords = ipw.Text(
            description="Vacancy coords:",
            style={"description_width": "150px"},
        )
        self.btn_defect_position.on_click(self._defect_position)
        self.btn_defect_position_vac.on_click(self._defect_position_vac)
        self.distortion_factor = ipw.BoundedFloatText(
            value=1.0,
            min=0.2,
            max=1.8,
            step=0.1,
            description="Distortion factor:",
            style={"description_width": "150px"},
        )
        self.btn_apply_bond_distortion = ipw.Button(
            description="Apply bond distortion",
            button_style="primary",
            disabled=False,
        )
        self.btn_apply_bond_distortion.on_click(self._apply_bond_distortion)
        self.selected_atoms = ipw.Text(
            description="Select atoms:", value="", style={"description_width": "150px"}
        )
        self.btn_selected_atoms = ipw.Button(
            description="From selection", style={"description_width": "150px"}
        )
        self.btn_selected_atoms.on_click(self._selected_atoms)
        self.wrong_syntax = ipw.HTML(
            value="""<i class="fa fa-times" style="color:red;font-size:2em;" ></i> wrong syntax""",
            layout={"visibility": "hidden"},
        )
        self.radial_cutoff = ipw.FloatText(
            description="Radial cutoff distance (Å):",
            value=3.0,
            style={"description_width": "150px"},
        )
        self.btn_apply_random_distorion_all = ipw.Button(
            description="Apply to all atoms",
            button_style="primary",
            style={"description_width": "250px"},
        )
        self.btn_apply_random_distorion_all.on_click(self._apply_random_distortion_all)
        self.btn_apply_random_distortion = ipw.Button(
            description="Apply locally",
            button_style="primary",
            style={"description_width": "250px"},
        )
        self.btn_apply_random_distortion.on_click(self._apply_random_distortion)
        super().__init__(
            children=[
                ipw.HTML(
                    """
                    <p>
                    Apply bond distortions or random displacements to explore metastable structures of point defects in solids.
                    This editor allows you to manipulate atomic positions and study defect behaviors in materials.
                    </p>
                    """
                ),
                ipw.HTML(
                    "<h4>Define defect position:</h4>",
                ),
                ipw.HTML(
                    """
                    <p>
                    To apply bond distortions, you can define either:<br>
                    - A <strong>defect atom</strong>: Select an atom in the viewer to mark it as a defect.<br>
                    - A <strong>vacancy position</strong>: Select the center of mass of the selection to define the vacancy position.
                    </p>
                    """
                ),
                ipw.HBox(
                    [
                        self.defect_position,
                        self.btn_defect_position,
                    ]
                ),
                ipw.HBox(
                    [
                        self.vacancy_coords,
                        self.btn_defect_position_vac,
                    ]
                ),
                ipw.HTML(
                    "<h4>Bond distortion around defect:</h4>",
                ),
                ipw.HTML(
                    """
                    <p>
                    Specify the number of neighboring atoms to apply the bond distortion.
                    Set a distortion factor to adjust atomic positions:<br>
                    - A value less than 1.0 shrinks the positions (e.g., 0.8 for a 20% reduction). <br>
                    - A value greater than 1.0 expands the positions (e.g., 1.8 for an 80% increase).
                    </p>
                    """
                ),
                ipw.VBox(
                    [
                        self.num_nearest_neighbours,
                        self.distortion_factor,
                        self.btn_apply_bond_distortion,
                    ]
                ),
                ipw.HTML(
                    "<h4>Random displacements to all atoms or selected ones:</h4>",
                ),
                ipw.HTML(
                    """
                    <p>
                    Define atoms for local random displacement around the defect or vacancy, and set the radial cutoff distance (in Å) to determine neighboring atoms for distance checks.
                    </p>
                    """
                ),
                ipw.HBox(
                    [
                        self.selected_atoms,
                        self.btn_selected_atoms,
                        self.wrong_syntax,
                    ]
                ),
                self.radial_cutoff,
                ipw.HBox(
                    [
                        self.btn_apply_random_distortion,
                        self.btn_apply_random_distorion_all,
                    ],
                ),
                self._status_message,
                ipw.HTML(
                    """
                    <p>
                    Uses the <a href="https://github.com/SMTG-Bham/ShakeNBreak" target="_blank">ShakeNBreak</a> package.
                    </p>
                    """
                ),
            ]
        )

    def sel2com(self):
        """Get center of mass of the selection."""
        if self.selection:
            com = np.average(
                self.structure[self.selection].get_scaled_positions(), axis=0
            )
        else:
            com = [0, 0, 0]

        return com

    def str2vec(self, string):
        return np.array(list(map(float, string.split())))

    def vec2str(self, vector):
        return (
            str(round(vector[0], 2))
            + " "
            + str(round(vector[1], 2))
            + " "
            + str(round(vector[2], 2))
        )

    def _defect_position_vac(self, _=None):
        """Define vaccancy coordinates."""
        self.vacancy_coords.value = self.vec2str(self.sel2com())

    def _defect_position(self, _=None):
        """Define vaccancy coordinates."""
        if len(self.selection) != 1:
            self._status_message.message = """
            <div class="alert alert-info">
            <strong>Please select one atom first.</strong>
            </div>
            """
        else:
            self.defect_position.value = list_to_string_range(self.selection)

    def _selected_atoms(self, _=None):
        """Selected atoms to displace."""
        self.selected_atoms.value = list_to_string_range(self.selection)

    def _apply_bond_distortion(self, _=None):
        if not self.defect_position.value and not self.vacancy_coords.value:
            self._status_message.message = """
            <div class="alert alert-info">
            <strong>Please select the position of the defect or vacancy.</strong>
            </div>
            """
        else:
            if self.defect_position.value == "":
                site_index = None
            else:
                site_index = int(self.defect_position.value)
            if self.vacancy_coords.value == "":
                frac_coords = None
            else:
                frac_coords = self.str2vec(self.vacancy_coords.value)

            self._apply_distortion(
                distort,
                num_nearest_neighbours=self.num_nearest_neighbours.value,
                site_index=site_index,
                distortion_factor=self.distortion_factor.value,
                frac_coords=frac_coords,
            )

    def _apply_random_distortion(self, _=None):
        if not self.defect_position.value and not self.vacancy_coords.value:
            self._status_message.message = """
            <div class="alert alert-info">
            <strong>Please select the position of the defect or vacancy.</strong>
            </div>
            """
        else:
            if self.defect_position.value == "":
                site_index = None
            else:
                site_index = int(self.defect_position.value)
            if self.vacancy_coords.value == "":
                frac_coords = None
            else:
                frac_coords = self.str2vec(self.vacancy_coords.value)

            active_atoms, syntax_ok = string_range_to_list(self.selected_atoms.value)
            if not active_atoms:
                active_atoms = None
            if not syntax_ok:
                self.wrong_syntax.layout.visibility = "visible"
            else:
                self.wrong_syntax.layout.visibility = "hidden"

                self._apply_distortion(
                    local_mc_rattle,
                    site_index=site_index,
                    frac_coords=frac_coords,
                    active_atoms=active_atoms,
                    nbr_cutoff=self.radial_cutoff.value,
                )

    def _apply_random_distortion_all(self, _=None):
        self._apply_distortion(
            rattle,
            active_atoms=None,
            nbr_cutoff=self.radial_cutoff.value,
        )

    def _apply_distortion(self, distortion_func, **kwargs):
        pymatgen_ase = AseAtomsAdaptor()
        pymatgen_structure = pymatgen_ase.get_structure(self.structure)
        struc_distorted = distortion_func(structure=pymatgen_structure, **kwargs)
        periodicity = self.structure.pbc
        if "distorted_structure" in struc_distorted:
            atoms = pymatgen_ase.get_atoms(struc_distorted["distorted_structure"])
        else:
            atoms = pymatgen_ase.get_atoms(struc_distorted)
        atoms.set_pbc(periodicity)
        self.structure = atoms
