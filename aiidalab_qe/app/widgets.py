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

# for AddTagsEditor
import ase
import ipywidgets as ipw
import numpy as np
import spglib
import traitlets
import traitlets as tl
from aiida import orm, plugins
from aiida.orm import CalcJobNode, load_node
from aiidalab_widgets_base import register_viewer_widget
from aiidalab_widgets_base.structures import _register_structure
from aiidalab_widgets_base.utils import (
    StatusHTML,
    list_to_string_range,
    string_range_to_list,
)
from IPython.display import HTML, Javascript, clear_output, display
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor

# defects
from shakenbreak.distortions import distort, local_mc_rattle

# trigger registration of the viewer widget:
from aiidalab_qe.app import node_view  # noqa: F401

StructureData = plugins.DataFactory("core.structure")


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
    structure_node = traitlets.Instance(orm.Data, allow_none=True, read_only=True)

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
        self.periodicity = ipw.RadioButtons(
            options=[
                "xyz",
                "xy",
                "x",
            ],
            value="xyz",
            description="Periodicty: ",
        )
        self.select_periodicity = ipw.Button(
            description="Select",
            button_style="primary",
            layout={"width": "initial"},
        )
        self.add_tags.on_click(self._add_tags)
        self.clear_tags.on_click(self._clear_tags)
        self.clear_all_tags.on_click(self._clear_all_tags)
        self.select_periodicity.on_click(self._select_periodicity)
        super().__init__(
            children=[
                ipw.HTML(
                    "<b>Adding a tag to atoms</b>",
                ),
                ipw.HBox([self.atom_selection, self.from_selection, self.tag]),
                ipw.HBox([self.add_tags, self.clear_tags, self.clear_all_tags]),
                ipw.HTML(
                    "<b>Define periodicity</b>",
                ),
                self.periodicity,
                self.select_periodicity,
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

    @_register_structure
    def _select_periodicity(self, _=None, atoms=None):
        """Select periodicity."""
        periodicity_options = {
            "xyz": (True, True, True),
            "xy": (True, True, False),
            "x": (True, False, False),
        }
        new_structure = deepcopy(self.structure)
        new_structure.pbc = periodicity_options[self.periodicity.value]
        self.structure = None
        self.structure = deepcopy(new_structure)


class DistorsionStructureEditor(ipw.VBox):
    structure = tl.Instance(ase.Atoms, allow_none=True)
    selection = tl.List(tl.Int)

    def __init__(self, title=""):
        self.title = title
        self._status_message = StatusHTML()
        self.num_nearest_neighbours = ipw.IntText(
            description="Num. neighbours atoms:",
            value=8,
            step=1,
            style={"description_width": "initial"},
            layout={"width": "initial"},
        )
        self.defect_position = ipw.Text(
            description="Defect atom:", layout={"width": "initial"}
        )
        # Only select one atom!
        self.btn_defect_position = ipw.Button(
            description="From selection",
            layout={"width": "initial"},
        )
        # Center of the selection
        self.btn_defect_position_vac = ipw.Button(
            description="From selection",
            layout={"width": "initial"},
        )
        self.vacancy_coords = ipw.Text(
            description="Vacancy coords:",
            layout={"width": "initial"},
            style={"description_width": "initial"},
        )
        self.btn_defect_position.on_click(self._defect_position)
        self.btn_defect_position_vac.on_click(self._defect_position_vac)
        self.distortion_factor = ipw.BoundedFloatText(
            value=1.0,
            min=0.2,
            max=1.8,
            step=0.1,
            description="Distortion Factor:",
            style={"description_width": "initial"},
            layout={"width": "initial"},
        )
        self.btn_apply_bond_distortion = ipw.Button(
            description="Apply Bond Distortion",
            button_style="primary",
            disabled=False,
            layout={"width": "initial"},
        )
        self.btn_apply_bond_distortion.on_click(self._apply_bond_distortion)
        self.selected_atoms = ipw.Text(
            description="Select atoms:",
            value="",
            style={"description_width": "initial"},
        )
        self.btn_selected_atoms = ipw.Button(
            description="From selection",
            layout={"width": "initial"},
        )
        self.btn_selected_atoms.on_click(self._selected_atoms)
        self.wrong_syntax = ipw.HTML(
            value="""<i class="fa fa-times" style="color:red;font-size:2em;" ></i> wrong syntax""",
            layout={"visibility": "hidden"},
        )
        self.radial_cutoff = ipw.FloatText(
            description="Radial cutoff distance (Å):",
            value=3,
            style={"description_width": "initial"},
            layout={"width": "initial"},
        )
        self.btn_apply_random_distortion = ipw.Button(
            description="Apply Random Displacement",
            button_style="primary",
            disabled=False,
            layout={"width": "initial"},
        )
        self.btn_apply_random_distortion.on_click(self._apply_random_distortion)
        super().__init__(
            children=[
                ipw.HTML(
                    "<b>Define defect position:</b>",
                ),
                ipw.HBox(
                    [
                        self.defect_position,
                        self.btn_defect_position,
                        self.vacancy_coords,
                        self.btn_defect_position_vac,
                    ]
                ),
                ipw.HTML(
                    "<b>Bond distortion around defect:</b>",
                ),
                ipw.HBox(
                    [
                        self.num_nearest_neighbours,
                        self.distortion_factor,
                    ]
                ),
                ipw.HTML(
                    "<b>Random displacements to all atoms or selected ones:</b>",
                ),
                ipw.HBox(
                    [
                        self.selected_atoms,
                        self.btn_selected_atoms,
                        self.wrong_syntax,
                        self.radial_cutoff,
                    ]
                ),
                ipw.HBox(
                    [
                        self.btn_apply_bond_distortion,
                        self.btn_apply_random_distortion,
                    ]
                ),
                self._status_message,
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

    @_register_structure
    def _apply_bond_distortion(self, _=None, atoms=None):
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
            pymatgen_ase = AseAtomsAdaptor()
            pymatgen_structure = pymatgen_ase.get_structure(atoms)
            struc_distorted = distort(
                structure=pymatgen_structure,
                num_nearest_neighbours=self.num_nearest_neighbours.value,
                site_index=site_index,
                distortion_factor=self.distortion_factor.value,
                frac_coords=frac_coords,
            )
            atoms = pymatgen_ase.get_atoms(struc_distorted["distorted_structure"])
            self.structure = atoms

    @_register_structure
    def _apply_random_distortion(self, _=None, atoms=None):
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

                pymatgen_ase = AseAtomsAdaptor()
                pymatgen_structure = pymatgen_ase.get_structure(atoms)
                struc_distorted = local_mc_rattle(
                    structure=pymatgen_structure,
                    site_index=site_index,
                    frac_coords=frac_coords,
                    active_atoms=active_atoms,
                    nbr_cutoff=self.radial_cutoff.value,
                )
                atoms = pymatgen_ase.get_atoms(struc_distorted)
                self.structure = atoms


class HubbardWidget(ipw.VBox):
    def __init__(self, input_structure=None):
        self.input_structure = input_structure
        self.hubbard = ipw.Checkbox(
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
                        self.hubbard,
                    ]
                ),
                self.hubbard_widget_out,
                self.eigen_values_widget_out,
            ]
        )
        self.hubbard.observe(self.toggle_hubbard_widgets, names="value")
        self.eigenvalues_label.observe(self.toggle_eigenvalues_widgets, names="value")

    def create_hubbard_widget(self):
        if self.input_structure is None:
            self.input_labels = []
        else:
            self.input_labels = self.input_structure.get_kind_names()
        widgets_list = []
        for label in self.input_labels:
            hbox_container = ipw.HBox()
            float_widget = ipw.BoundedFloatText(
                description=label, min=0, max=10, step=0.1, value=0.0
            )
            hbox_container.children = [float_widget]
            widgets_list.append(hbox_container)
        hubbard_widget = ipw.VBox(
            [ipw.HTML("Define U value")] + widgets_list + [self.eigenvalues_label]
        )
        return hubbard_widget

    def create_eigenvalues_widget(self):
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
                num_states = 5
            if kind[2].is_lanthanoid:
                num_states = 7
            if kind[2].is_transition_metal or kind[2].is_lanthanoid:
                widgets_list_up = []
                widgets_list_down = []
                for i in range(num_states):
                    eigenvalues_up = ipw.Dropdown(
                        description=f"{i+1}",
                        options=["-1", "0", "1"],
                        layout=ipw.Layout(width="150px"),
                    )
                    eigenvalues_down = ipw.Dropdown(
                        description=f"{i+1}",
                        options=["-1", "0", "1"],
                        layout=ipw.Layout(width="150px"),
                    )
                    widgets_list_up.append(eigenvalues_up)
                    widgets_list_down.append(eigenvalues_down)

                row_up = ipw.HBox(
                    children=[
                        ipw.Label(
                            "Up:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="40px"
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
                                justify_content="flex-start", width="40px"
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
                                justify_content="flex-start", width="40px"
                            ),
                        ),
                        eigenvalues_container,
                    ]
                )
                kind_list.append(kind_container)
        occup_kinds_widget = ipw.VBox(kind_list)

        return occup_kinds_widget

    def update_widgets(self, change):
        self.input_structure = change["new"]
        self.input_labels = self.input_structure.get_kind_names()
        self.hubbard_widget = self.create_hubbard_widget()
        if self.hubbard.value:
            with self.hubbard_widget_out:
                clear_output()
                display(self.hubbard_widget)
        self.eigen_values_widget = self.create_eigenvalues_widget()
        if self.eigenvalues_label.value:
            with self.eigen_values_widget_out:
                clear_output()
                display(self.eigen_values_widget)

    # def update_hubbard_widgets(self, change):
    #     self.input_labels = self.input_structure.get_kind_names()
    #     self.hubbard_widget = self.create_hubbard_widget()
    #     if self.hubbard.value:
    #         with self.hubbard_widget_out:
    #             clear_output()
    #             display(self.hubbard_widget)

    def toggle_hubbard_widgets(self, change):
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
        if change["new"]:
            with self.eigen_values_widget_out:
                clear_output()
                display(self.eigen_values_widget)
        else:
            with self.eigen_values_widget_out:
                clear_output()

    def _get_hubbard_u(self):
        hubbard_u = {}
        for index, label in enumerate(self.input_labels):
            value_hubbard = self.hubbard_widget.children[index + 1].children[0].value
            if value_hubbard != 0:
                hubbard_u[label] = (
                    self.hubbard_widget.children[index + 1].children[0].value
                )
        return hubbard_u

    def _get_starting_ns_eigenvalue(self) -> list:
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
                            # starting_ns_eigenvalue.append([j+1, i+1, kind[0], value_eigenvalue])
                            starting_ns_eigenvalue.append(
                                [j + 1, i + 1, kind[1], value_eigenvalue]
                            )

        return starting_ns_eigenvalue

    @property
    def hubbard_dict(self) -> dict:
        if self.hubbard.value:
            hubbard_dict = {
                "hubbard_u": self._get_hubbard_u(),
                "lda_plus_u": True,
                "U_projection_type": "ortho-atomic",
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


class BasicCellEditor(ipw.VBox):
    """Widget that allows for the basic cell editing."""

    structure = tl.Instance(ase.Atoms, allow_none=True)

    def __init__(self, title="Cell transform"):
        self.title = title
        self._status_message = StatusHTML()

        # cell transfor opration widget
        primitive_cell = ipw.Button(
            description="Convert to primitive cell",
            layout={"width": "initial"},
        )
        primitive_cell.on_click(self._to_primitive_cell)

        conventional_cell = ipw.Button(
            description="Convert to conventional cell",
            layout={"width": "initial"},
        )
        conventional_cell.on_click(self._to_conventional_cell)
        # change PR 372
        cell_parameters = (
            self.structure.get_cell_lengths_and_angles()
            if self.structure
            else [1, 0, 0, 0, 0, 0]
        )
        self.cell_parameters = ipw.HBox(
            [
                ipw.VBox(
                    [
                        ipw.HTML(
                            description=["a(Å)", "b(Å)", "c(Å)", "α", "β", "γ"][i],
                            layout={"width": "30px"},
                        ),
                        ipw.FloatText(
                            value=cell_parameters[i], layout={"width": "100px"}
                        ),
                    ]
                )
                for i in range(6)
            ]
        )
        #
        self.cell_transformation = ipw.VBox(
            [
                ipw.HBox(
                    [
                        ipw.IntText(value=1 if i == j else 0, layout={"width": "60px"})
                        for j in range(3)
                    ]
                    + [ipw.FloatText(value=0, layout={"width": "60px"})]
                )
                for i in range(3)
            ]
        )
        apply_cell_parameters = ipw.Button(description="Apply cell parameters")
        apply_cell_parameters.on_click(self.apply_cell_parameters)
        apply_cell_transformation = ipw.Button(description="Apply transformation")
        apply_cell_transformation.on_click(self.apply_cell_transformation)
        # change ends
        super().__init__(
            children=[
                ipw.HBox(
                    [
                        primitive_cell,
                        conventional_cell,
                    ],
                ),
                self._status_message,
                ipw.VBox(  # Change here
                    [
                        ipw.HTML(
                            "<b>Cell parameters:</b>",
                            layout={"margin": "20px 0px 10px 0px"},
                        ),
                        self.cell_parameters,
                        apply_cell_parameters,
                    ],
                    layout={"margin": "0px 0px 0px 20px"},
                ),
                ipw.VBox(
                    [
                        ipw.HTML(
                            "<b>Cell Transformation:</b>",
                            layout={"margin": "20px 0px 10px 0px"},
                        ),
                        self.cell_transformation,
                        apply_cell_transformation,
                    ],
                    layout={"margin": "0px 0px 0px 20px"},
                ),  # change here
            ],
        )

    @_register_structure
    def _to_primitive_cell(self, _=None, atoms=None):
        atoms = self._to_standard_cell(atoms, to_primitive=True)

        self.structure = atoms

    @_register_structure
    def _to_conventional_cell(self, _=None, atoms=None):
        atoms = self._to_standard_cell(atoms, to_primitive=False)

        self.structure = atoms

    @staticmethod
    def _to_standard_cell(
        structure: ase.Atoms, to_primitive=False, no_idealize=False, symprec=1e-5
    ):
        """The `standardize_cell` method from spglib and apply to ase.Atoms"""
        lattice = structure.get_cell()
        positions = structure.get_scaled_positions()
        numbers = structure.get_atomic_numbers()

        cell = (lattice, positions, numbers)

        # operation
        lattice, positions, numbers = spglib.standardize_cell(
            cell, to_primitive=to_primitive, no_idealize=no_idealize, symprec=symprec
        )

        return ase.Atoms(
            cell=lattice,
            scaled_positions=positions,
            numbers=numbers,
            pbc=[True, True, True],
        )

    @tl.observe("structure")  # change here
    def _observe_structure(self, change):
        """Update cell after the structure has been modified."""
        if change["new"] is not None:
            cell_parameters = change["new"].get_cell_lengths_and_angles()
            for i in range(6):
                self.cell_parameters.children[i].children[1].value = round(
                    cell_parameters[i], 4
                )
        else:
            for i in range(6):
                self.cell_parameters.children[i].children[1].value = 0

    @_register_structure
    def apply_cell_parameters(self, _=None, atoms=None):
        from ase.cell import Cell

        # only update structure when atoms is not None.
        cell_parameters = [
            self.cell_parameters.children[i].children[1].value for i in range(6)
        ]
        if atoms is not None:
            atoms.cell = Cell.fromcellpar(cell_parameters)
            self.structure = atoms

    @_register_structure
    def apply_cell_transformation(self, _=None, atoms=None):
        from ase.build import make_supercell

        # only update structure when atoms is not None.
        if atoms is not None:
            mat = np.zeros((3, 3))
            translate = np.zeros(3)
            for i in range(3):
                translate[i] = self.cell_transformation.children[i].children[3].value
                for j in range(3):
                    mat[i][j] = self.cell_transformation.children[i].children[j].value
            # transformation matrix may not work due to singularity, or
            # the number of generated atoms is not correct
            try:
                atoms = make_supercell(atoms, mat)
            except Exception as e:
                self._status_message.message = """
            <div class="alert alert-info">
            <strong>The transformation matrix is wrong! {}</strong>
            </div>
            """.format(
                    e
                )
                return
            # translate
            atoms.translate(-atoms.cell.array.dot(translate))
            self.structure = atoms  # change here
