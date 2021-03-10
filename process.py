"""Widgets for the monitoring of processes."""
import os
from dataclasses import dataclass
from threading import Thread, Event, Lock

import traitlets
import ipywidgets as ipw
from aiida.cmdline.utils.query.calculation import CalculationQueryBuilder
from aiida.orm import ProcessNode
from aiida.orm import load_node


def get_calc_job_output(process):
    from aiidalab_widgets_base.process import get_running_calcs

    previous_calc_id = None
    num_lines = 0

    while not process.is_sealed:
        calc = None
        for calc in get_running_calcs(process):
            if calc.id == previous_calc_id:
                break
        else:
            if calc:
                previous_calc_id = calc.id

        if calc and "remote_folder" in calc.outputs:
            f_path = os.path.join(
                calc.outputs.remote_folder.get_remote_path(),
                calc.attributes["output_filename"],
            )
            if os.path.exists(f_path):
                with open(f_path) as fobj:
                    new_lines = fobj.readlines()[num_lines:]
                    num_lines += len(new_lines)
                    yield from new_lines


class ProgressBarWidget(ipw.VBox):
    """A bar showing the proggress of a process."""

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        self.correspondance = {
            None: (0, "warning"),
            "created": (0, "info"),
            "running": (1, "info"),
            "waiting": (1, "info"),
            "killed": (2, "danger"),
            "excepted": (2, "danger"),
            "finished": (2, "success"),
        }
        self.bar = ipw.IntProgress(  # pylint: disable=blacklisted-name
            value=0,
            min=0,
            max=2,
            step=1,
            bar_style="warning",  # 'success', 'info', 'warning', 'danger' or ''
            orientation="horizontal",
            layout=ipw.Layout(width="auto"),
        )
        self.state = ipw.HTML(
            description="Calculation state:",
            value="",
            style={"description_width": "100px"},
        )
        super().__init__(children=[self.state, self.bar], **kwargs)

    @traitlets.observe("process")
    def update(self, _=None):
        """Update the bar."""
        self.bar.value, self.bar.bar_style = self.correspondance[self.current_state]
        if self.current_state is None:
            self.state.value = "N/A"
        else:
            self.state.value = self.current_state.capitalize()

    @property
    def current_state(self):
        if self.process is not None:
            return self.process.process_state.value


class WorkChainSelector(ipw.HBox):

    # The PK of a 'aiida.workflows:quantumespresso.pw.bands' WorkChainNode.
    value = traitlets.Int(allow_none=True)

    # When this trait is set to a positive value, the work chains are automatically
    # refreshed every `auto_refresh_interval` seconds.
    auto_refresh_interval = traitlets.Int()  # seconds

    # Indicate whether the widget is currently updating the work chain options.
    busy = traitlets.Bool(read_only=True)

    # Note: We use this class as a singleton to reset the work chains selector
    # widget to its default stage (no work chain selected), because we cannot
    # use `None` as setting the widget's value to None will lead to "no selection".
    _NO_PROCESS = object()

    FMT_WORKCHAIN = "{wc.pk:6}{wc.ctime:>10}\t{wc.state:<16}\t{wc.formula}"

    def __init__(self, **kwargs):
        self.work_chains_selector = ipw.Dropdown(
            description="WorkChain",
            options=[("New calculation...", self._NO_PROCESS)],
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
        )
        ipw.dlink(
            (self.work_chains_selector, "value"),
            (self, "value"),
            transform=lambda pk: None if pk is self._NO_PROCESS else pk,
        )

        self.refresh_work_chains_button = ipw.Button(description="Refresh")
        self.refresh_work_chains_button.on_click(self.refresh_work_chains)

        self._refresh_lock = Lock()
        self._refresh_thread = None
        self._stop_refresh_thread = Event()
        self._update_auto_refresh_thread_state()

        super().__init__(
            children=[self.work_chains_selector, self.refresh_work_chains_button],
            **kwargs,
        )

    @dataclass
    class WorkChainData:
        pk: int
        ctime: str
        state: str
        formula: str

    @classmethod
    def find_work_chains(cls):
        builder = CalculationQueryBuilder()
        filters = builder.get_filters(
            process_label="PwBandsWorkChain",
        )
        query_set = builder.get_query_set(
            filters=filters,
            order_by={"ctime": "desc"},
        )
        projected = builder.get_projected(
            query_set, projections=["pk", "ctime", "state"]
        )

        for process in projected[1:]:
            pk = process[0]
            formula = load_node(pk).inputs.structure.get_formula()
            yield cls.WorkChainData(formula=formula, *process)

    @traitlets.default("busy")
    def _default_busy(self):
        return True

    @traitlets.observe("busy")
    def _observe_busy(self, change):
        for child in self.children:
            child.disabled = change["new"]

    def refresh_work_chains(self, _=None):
        with self._refresh_lock:
            try:
                self.set_trait("busy", True)  # disables the widget

                with self.hold_trait_notifications():
                    # We need to restore the original value, because it may be reset due to this issue:
                    # https://github.com/jupyter-widgets/ipywidgets/issues/2230
                    original_value = self.work_chains_selector.value

                    self.work_chains_selector.options = [
                        ("New calculation...", self._NO_PROCESS)
                    ] + [
                        (self.FMT_WORKCHAIN.format(wc=wc), wc.pk)
                        for wc in self.find_work_chains()
                    ]

                    self.work_chains_selector.value = original_value
            finally:
                self.set_trait("busy", False)  # reenable the widget

    def _auto_refresh_loop(self):
        while True:
            self.refresh_work_chains()
            if self._stop_refresh_thread.wait(timeout=self.auto_refresh_interval):
                break

    def _update_auto_refresh_thread_state(self):
        if self.auto_refresh_interval > 0 and self._refresh_thread is None:
            # start thread
            self._stop_refresh_thread.clear()
            self._refresh_thread = Thread(target=self._auto_refresh_loop)
            self._refresh_thread.start()

        elif self.auto_refresh_interval <= 0 and self._refresh_thread is not None:
            # stop thread
            self._stop_refresh_thread.set()
            self._refresh_thread.join(timeout=30)
            self._refresh_thread = None

    @traitlets.default("auto_refresh_interval")
    def _default_auto_refresh_interval(self):
        return 10  # seconds

    @traitlets.observe("auto_refresh_interval")
    def _observe_auto_refresh_interval(self, change):
        if change["new"] != change["old"]:
            self._update_auto_refresh_thread_state()

    @traitlets.observe("value")
    def _observe_value(self, change):
        if change["old"] == change["new"]:
            return

        new = self._NO_PROCESS if change["new"] is None else change["new"]

        if new not in {pk for _, pk in self.work_chains_selector.options}:
            self.refresh_work_chains()

        self.work_chains_selector.value = new
