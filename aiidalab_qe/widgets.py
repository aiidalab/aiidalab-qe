"""Widgets for the QE app.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""

import base64
from tempfile import NamedTemporaryFile
from threading import Event
from threading import Lock
from threading import Thread
from queue import Queue

import ipywidgets as ipw
import traitlets
from aiida.orm import CalcJobNode
from aiida.orm import Node
from aiidalab_widgets_base import viewer
from aiidalab_widgets_base import register_viewer_widget
from IPython.display import clear_output
from IPython.display import display

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


class ResourceSelectionWidget(ipw.HBox):
    """Widget for the selection of compute (CPU) resources."""

    resource_selection_title = ipw.HTML("""<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Resources</h4>
    </div>""")
    resource_selection_prompt = ipw.HTML(
        "Select the compute resources for this calculation."
    )
    resource_selection_help = ipw.HTML(
        """<div style="line-height:120%; padding-top:0px">
        <p style="padding-bottom:10px">
        In this tab you can select the amount of resources you want to use for the calculation.
        There is no general rule of thumb on how to select the appropriate number of nodes and cores. In general:
        </p>
        <ul>
        <li>Increase the number of nodes if you run out of memory for larger structures.</li>
        <li>Increase the number of nodes and cores if you want to reduce the total runtime.</li>
        </ul>
        <p>However, specifying the optimal configuration of resources is a complex issue and
        simply increasing either cores or nodes may not have the desired effect.</p></div>"""
    )

    def __init__(self, **kwargs):
        extra = {
            "style": {"description_width": "150px"},
            "layout": {"max_width": "200px"},
        }
        self.number_of_nodes = ipw.BoundedIntText(
            value=1, step=1, min=1, description="# nodes", disabled=False, **extra
        )
        self.cpus_per_node = ipw.BoundedIntText(
            value=1, step=1, min=1, description="# cpus per node", **extra
        )
        self.total_num_cpus = ipw.BoundedIntText(
            value=1, step=1, min=1, description="# total cpus", disabled=True, **extra
        )

        # Update the total # of CPUs int text:
        self.number_of_nodes.observe(self._update_total_num_cpus, "value")
        self.cpus_per_node.observe(self._update_total_num_cpus, "value")

        super().__init__(
            children=[
                ipw.VBox(
                    children=[
                        self.resource_selection_title,
                        ipw.HBox(
                            children=[
                                self.resource_selection_help,
                                ipw.VBox(
                                    children=[
                                        self.number_of_nodes,
                                        self.cpus_per_node,
                                        self.total_num_cpus,
                                    ],
                                    layout=ipw.Layout(min_width="310px"),
                                ),
                            ]
                        )
                    ]
                )
            ]
        )

    def _update_total_num_cpus(self, change):
        self.total_num_cpus.value = (
            self.number_of_nodes.value * self.cpus_per_node.value
        )

    def reset(self):
        with self.hold_trait_notifications():
            self.number_of_nodes.value = 1
            self.cpus_per_node.value = 1
