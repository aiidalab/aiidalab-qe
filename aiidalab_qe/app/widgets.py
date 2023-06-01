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
from aiida.orm import CalcJobNode, load_node
from aiidalab_widgets_base import register_viewer_widget
from aiidalab_widgets_base.utils import (
    StatusHTML,
    list_to_string_range,
    string_range_to_list,
)
from IPython.display import HTML, Javascript, display

# trigger registration of the viewer widget:
from aiidalab_qe.app import node_view  # noqa: F401

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
        self._html = ipw.HTML()
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


@register_viewer_widget("process.calculation.calcjob.CalcJobNode.")
class CalcJobNodeViewerWidget(ipw.VBox):
    def __init__(self, calcjob, **kwargs):
        self.calcjob = calcjob
        self.output_follower = CalcJobOutputFollower()
        self.log_output = LogOutputWidget()

        self.output_follower.calcjob_uuid = self.calcjob.uuid
        self.output_follower.observe(self._observe_output_follower_lineno, ["lineno"])

        super().__init__(
            [ipw.HTML(f"CalcJob: {self.calcjob}"), self.log_output], **kwargs
        )

    def _observe_output_follower_lineno(self, _):
        with self.hold_trait_notifications():
            self.log_output.filename = self.output_follower.filename
            self.log_output.value = "\n".join(self.output_follower.output)


class ResourceSelectionWidget(ipw.VBox):
    """Widget for the selection of compute resources."""

    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Resources</h4>
    </div>"""
    )
    prompt = ipw.HTML(
        """<div style="line-height:120%; padding-top:0px">
        <p style="padding-bottom:10px">
        Specify the resources to use for the pw.x calculation.
        </p></div>"""
    )

    def __init__(self, **kwargs):
        extra = {
            "style": {"description_width": "150px"},
            "layout": {"min_width": "180px"},
        }
        self.num_nodes = ipw.BoundedIntText(
            value=1, step=1, min=1, max=1000, description="Nodes", **extra
        )
        self.num_cpus = ipw.BoundedIntText(
            value=1, step=1, min=1, description="CPUs", **extra
        )

        super().__init__(
            children=[
                self.title,
                ipw.HBox(
                    children=[self.prompt, self.num_nodes, self.num_cpus],
                    layout=ipw.Layout(justify_content="space-between"),
                ),
            ]
        )

    def reset(self):
        self.num_nodes.value = 1
        self.num_cpus.value = 1


class ParallelizationSettings(ipw.VBox):
    """Widget for setting the parallelization settings."""

    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Parallelization</h4>
    </div>"""
    )
    prompt = ipw.HTML(
        """<div style="line-height:120%; padding-top:0px">
        <p style="padding-bottom:10px">
        Specify the number of k-points pools for the calculations.
        </p></div>"""
    )

    def __init__(self, **kwargs):
        extra = {
            "style": {"description_width": "150px"},
            "layout": {"min_width": "180px"},
        }
        self.npools = ipw.BoundedIntText(
            value=1, step=1, min=1, max=128, description="Number of k-pools", **extra
        )
        super().__init__(
            children=[
                self.title,
                ipw.HBox(
                    children=[self.prompt, self.npools],
                    layout=ipw.Layout(justify_content="space-between"),
                ),
            ]
        )

    def reset(self):
        self.npools.value = 1


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

    def __init__(self, title=""):
        self.title = title
        self._status_message = StatusHTML()
        self.atom_selection = ipw.Text(
            description="Define kind", value="", layout={"width": "initial"}
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

        self.clear_tags = ipw.Button(
            description="Clear tags",
            button_style="primary",
            layout={"width": "initial"},
        )
        self.clear_all_tags = ipw.Button(
            description="Clear all tags",
            button_style="primary",
            layout={"width": "initial"},
        )
        self.add_tags.on_click(self._add_tags)
        self.clear_tags.on_click(self._clear_tags)
        self.clear_all_tags.on_click(self._clear_all_tags)
        super().__init__(
            children=[
                ipw.HTML(
                    "<b>Adding a tag to atoms</b>",
                ),
                ipw.HBox([self.atom_selection, self.from_selection, self.tag]),
                ipw.HBox([self.add_tags, self.clear_tags, self.clear_all_tags]),
                self._status_message,
            ]
        )

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
            selection = string_range_to_list(self.atom_selection.value)
            new_structure = deepcopy(self.structure)
            if new_structure.get_tags() == []:
                new_tags = np.zeros(len(new_structure))
            else:
                new_tags = new_structure.get_tags()
            new_tags[selection] = self.tag.value
            new_structure.set_tags(new_tags)
            self.structure = None
            self.structure = deepcopy(new_structure)
            self.input_selection = None
            self.input_selection = deepcopy(self.selection)

    def _clear_tags(self, _=None):
        """Clear tags from selected atoms."""
        if not self.atom_selection.value:
            self._status_message.message = """
            <div class="alert alert-info">
            <strong>Please select atoms first.</strong>
            </div>
            """
        else:
            selection = string_range_to_list(self.atom_selection.value)
            new_structure = deepcopy(self.structure)
            new_tags = new_structure.get_tags()
            new_tags[selection] = 0
            new_structure.set_tags(new_tags)
            self.structure = None
            self.structure = deepcopy(new_structure)
            self.input_selection = None
            self.input_selection = deepcopy(self.selection)

    def _clear_all_tags(self, _=None):
        """Clear all tags."""
        new_structure = deepcopy(self.structure)
        new_tags = np.zeros(len(new_structure))
        new_structure.set_tags(new_tags)
        self.structure = None
        self.structure = deepcopy(new_structure)
        self.input_selection = None
        self.input_selection = deepcopy(self.selection)
