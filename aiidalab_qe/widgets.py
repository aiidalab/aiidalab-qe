"""Widgets for the QE app.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""

import base64
from queue import Queue
from tempfile import NamedTemporaryFile
from threading import Event, Lock, Thread
from time import time

import ipywidgets as ipw
import traitlets
from aiida.orm import CalcJobNode, Node
from aiidalab_widgets_base import register_viewer_widget, viewer
from IPython.display import clear_output, display

# trigger registration of the viewer widget:
from aiidalab_qe import node_view  # noqa: F401

__all__ = [
    "CalcJobOutputFollower",
    "LogOutputWidget",
    "NodeViewWidget",
]


class LogOutputWidget(ipw.VBox):

    value = traitlets.Tuple(traitlets.Unicode(), traitlets.Unicode())

    def __init__(self, num_min_lines=10, max_output_height="200px", **kwargs):
        self._num_min_lines = num_min_lines

        self._filename = ipw.Text(
            description="Filename:",
            placeholder="No file selected.",
            disabled=True,
            layout=ipw.Layout(flex="1 1 auto", width="auto"),
        )

        self._output = ipw.HTML(layout=ipw.Layout(min_width="50em"))
        self._output_container = ipw.VBox(
            children=[self._output], layout=ipw.Layout(max_height=max_output_height)
        )
        self._download_link = ipw.HTML(layout=ipw.Layout(width="auto"))

        self._refresh_output()
        super().__init__(
            children=[
                self._output_container,
                ipw.HBox(
                    children=[self._filename, self._download_link],
                    layout=ipw.Layout(min_height="25px"),
                ),
            ],
            **kwargs,
        )

    @traitlets.default("value")
    def _default_value(self):
        if self._num_min_lines > 0:
            return "", "\n" * self._num_min_lines

    @traitlets.observe("value")
    def _refresh_output(self, _=None):
        filename, loglines = self.value
        with self.hold_trait_notifications():
            self._filename.value = filename
            self._output.value = self._format_output(loglines)

            payload = base64.b64encode(loglines.encode()).decode()
            html_download = f'<a download="{filename}" href="data:text/plain;base64,{payload}" target="_blank">Download</a>'
            self._download_link.value = html_download

    style = "background-color: #253239; color: #cdd3df; line-height: normal"

    def _format_output(self, text):
        lines = text.splitlines()

        # Add empty lines to reach the minimum number of lines.
        lines += [""] * max(0, self._num_min_lines - len(lines))

        # Replace empty lines with single white space to ensure that they are actually shown.
        lines = [line if len(line) > 0 else " " for line in lines]

        # Replace the first line if there is no output whatsoever
        if len(text.strip()) == 0 and len(lines) > 0:
            lines[0] = "[waiting for output]"

        text = "\n".join(lines)
        return f"""<pre style="{self.style}">{text}</pre>"""


class CalcJobOutputFollower(traitlets.HasTraits):

    calcjob = traitlets.Instance(CalcJobNode, allow_none=True)
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

    @traitlets.observe("calcjob")
    def _observe_calcjob(self, change):
        try:
            if change["old"].pk == change["new"].pk:
                # Old and new process are identical.
                return
        except AttributeError:
            pass

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
                    target=self._follow_output, args=(change["new"],)
                )
                self._follow_output_thread.start()

    def _follow_output(self, calcjob):
        """Monitor calcjob and orchestrate pushing and pulling of output."""
        self._pull_thread = Thread(target=self._pull_output, args=(calcjob,))
        self._pull_thread.start()
        self._push_thread = Thread(target=self._push_output, args=(calcjob,))
        self._push_thread.start()

    def _fetch_output(self, calcjob):
        assert isinstance(calcjob, CalcJobNode)
        if "remote_folder" in calcjob.outputs:
            try:
                fn_out = calcjob.attributes["output_filename"]
                self.filename = fn_out
                with NamedTemporaryFile() as tmpfile:
                    calcjob.outputs.remote_folder.getfile(fn_out, tmpfile.name)
                    return tmpfile.read().decode().splitlines()
            except OSError:
                return list()
        else:
            return list()

    _EOF = None

    def _push_output(self, calcjob, delay=0.2):
        """Push new log lines onto the queue."""
        lineno = 0
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
                    break

    def _pull_output(self, calcjob):
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

        self.output_follower.observe(self._observe_output_follower_lineno, ["lineno"])
        self.output_follower.calcjob = self.calcjob

        super().__init__(
            [ipw.HTML(f"CalcJob: {self.calcjob}"), self.log_output], **kwargs
        )

    def _observe_output_follower_lineno(self, change):
        restrict_num_lines = None if self.calcjob.is_sealed else -10
        new_lines = "\n".join(self.output_follower.output[restrict_num_lines:])
        self.log_output.value = self.output_follower.filename, new_lines


class NodeViewWidget(ipw.VBox):

    node = traitlets.Instance(Node, allow_none=True)

    def __init__(self, **kwargs):
        self._output = ipw.Output()
        super().__init__(children=[self._output], **kwargs)

    @traitlets.observe("node")
    def _observe_node(self, change):
        if change["new"] != change["old"]:
            with self._output:
                clear_output()
                if change["new"]:
                    display(viewer(change["new"]))


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
        Specify the number of MPI tasks for this calculation.
        In general, larger structures will require a larger number of tasks.
        </p></div>"""
    )

    def __init__(self, **kwargs):
        extra = {
            "style": {"description_width": "150px"},
            # "layout": {"max_width": "200px"},
            "layout": {"min_width": "310px"},
        }

        self.num_mpi_tasks = ipw.BoundedIntText(
            value=1, step=1, min=1, description="# MPI tasks", **extra
        )

        super().__init__(
            children=[
                self.title,
                ipw.HBox(children=[self.prompt, self.num_mpi_tasks]),
            ]
        )

    def reset(self):
        self.num_mpi_tasks.value = 1


class ProgressBar(ipw.HBox):
    class AnimationRate(float):
        def __eq__(self, other):
            return type(self) == type(other) and super().__eq__(other)

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
