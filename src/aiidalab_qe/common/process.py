"""Widgets related to process management."""

from dataclasses import make_dataclass

import ipywidgets as ipw
import traitlets as tl


class WorkChainSelector(ipw.HBox):
    """A widget to select a WorkChainNode of a given process label.

    To use it, subclass it and set the `process_label` attribute to the desired.

    If you want to display additional information about the work chain, set the
    `extra_fields` attribute to a list of tuples, where each tuple contains the
    name of the field and the type of the field. The field names must match the
    names of the output keys of the `parse_extra_info` method.
    """

    # The PK of a WorkChainNode.
    value = tl.Int(allow_none=True)

    # Indicate whether the widget is currently updating the work chain options.
    busy = tl.Bool(read_only=True)

    # Note: We use this class as a singleton to reset the work chains selector
    # widget to its default stage (no work chain selected), because we cannot
    # use `None` as setting the widget's value to None will lead to "no selection".
    _NO_PROCESS = object()

    BASE_FMT_WORKCHAIN = "{wc.pk:6}{wc.ctime:>10}\t{wc.state:<16}\t{wc.label:<16}"

    projections = ["pk", "ctime", "state", "label"]
    BASE_FIELDS = [("pk", int), ("ctime", str), ("state", str), ("label", str)]
    extra_fields = None

    def __init__(self, process_label, **kwargs):
        self.process_label = process_label
        self.work_chains_prompt = ipw.HTML(
            "<b>Select computed workflow or start a new one:</b>&nbsp;"
        )
        self.work_chains_selector = ipw.Dropdown(
            options=[("New workflow...", self._NO_PROCESS)],
            layout=ipw.Layout(min_width="360px", flex="1 1 auto"),
        )
        ipw.dlink(
            (self.work_chains_selector, "value"),
            (self, "value"),
            transform=lambda pk: None if pk is self._NO_PROCESS else pk,
        )

        if self.extra_fields is not None:
            fmt_extra = "".join([f"{{wc.{field[0]}}}" for field in self.extra_fields])
            self.fmt_workchain = self.BASE_FMT_WORKCHAIN + "\t" + fmt_extra
        else:
            self.fmt_workchain = self.BASE_FMT_WORKCHAIN

        self.new_work_chains_button = ipw.Button(
            description="New",
            tooltip="Start a new computation",
            button_style="success",
            icon="plus-circle",
            layout=ipw.Layout(width="auto"),
        )
        self.new_work_chains_button.on_click(self._on_click_new_work_chain)

        self.refresh_work_chains_button = ipw.Button(
            description="Refresh",
            tooltip="Refresh the list of computed workflows",
            button_style="primary",
            icon="refresh",
            layout=ipw.Layout(width="auto"),
        )
        self.refresh_work_chains_button.on_click(self.refresh_work_chains)

        super().__init__(
            children=[
                self.work_chains_prompt,
                self.work_chains_selector,
                self.new_work_chains_button,
                self.refresh_work_chains_button,
            ],
            **kwargs,
        )

        self.refresh_work_chains()
        # the following is needed to disable the button.

    def parse_extra_info(self, _pk: int) -> dict:
        """Parse extra information about the work chain."""
        return {}

    def find_work_chains(self):
        from aiida.tools.query.calculation import CalculationQueryBuilder

        builder = CalculationQueryBuilder()
        filters = builder.get_filters(
            process_label=self.process_label,
        )
        query_set = builder.get_query_set(
            filters=filters,
            order_by={"ctime": "desc"},
        )
        projected = builder.get_projected(
            query_set,
            projections=self.projections,
        )

        for result in projected[1:]:
            process_info = dict(zip(self.projections, result))

            if self.extra_fields is not None:
                pk = process_info["pk"]
                extra_info = self.parse_extra_info(pk)

                yield make_dataclass("WorkChain", self.BASE_FIELDS + self.extra_fields)(
                    **process_info, **extra_info
                )
            else:
                yield make_dataclass("WorkChain", self.BASE_FIELDS)(**process_info)

    @tl.default("busy")
    def _default_busy(self):
        return True

    @tl.observe("busy")
    def _observe_busy(self, change):
        for child in self.children:
            child.disabled = change["new"]

    def refresh_work_chains(self, _=None):
        try:
            self.set_trait("busy", True)  # disables the widget

            with self.hold_trait_notifications():
                # We need to restore the original value, because it may be reset due to this issue:
                # https://github.com/jupyter-widgets/ipywidgets/issues/2230
                # It is fixed in ipywidgets 8.0.0, but we still support ipywidgets 7.x.
                original_value = self.work_chains_selector.value

                self.work_chains_selector.options = [
                    ("New workflow...", self._NO_PROCESS)
                ] + [
                    (self.fmt_workchain.format(wc=wc), wc.pk)
                    for wc in self.find_work_chains()
                ]

                self.work_chains_selector.value = original_value
        finally:
            self.set_trait("busy", False)  # reenable the widget

    def _on_click_new_work_chain(self, _=None):
        self.refresh_work_chains()
        self.work_chains_selector.value = self._NO_PROCESS

    @tl.observe("value")
    def _observe_value(self, change):
        if change["old"] == change["new"]:
            return

        new = self._NO_PROCESS if change["new"] is None else change["new"]

        if new not in {pk for _, pk in self.work_chains_selector.options}:
            self.refresh_work_chains()

        self.work_chains_selector.value = new


class QeAppWorkChainSelector(WorkChainSelector):
    # extra_fields = [("formula", str), ("relax_info", str), ("properties_info", str)]

    def __init__(self, **kwargs):
        super().__init__(process_label="QeAppWorkChain", **kwargs)
