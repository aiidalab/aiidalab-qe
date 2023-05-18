"""Widgets related to process management."""
from dataclasses import make_dataclass

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiida.tools.query.calculation import CalculationQueryBuilder


class WorkChainSelector(ipw.HBox):
    # The PK of a 'aiida.workflows:quantumespresso.pw.bands' WorkChainNode.
    value = tl.Int(allow_none=True)

    # Indicate whether the widget is currently updating the work chain options.
    busy = tl.Bool(read_only=True)

    # Note: We use this class as a singleton to reset the work chains selector
    # widget to its default stage (no work chain selected), because we cannot
    # use `None` as setting the widget's value to None will lead to "no selection".
    _NO_PROCESS = object()

    FMT_WORKCHAIN = "{wc.pk:6}{wc.ctime:>10}\t{wc.state:<16}"

    BASE_FIELDS = [("pk", int), ("ctime", str), ("state", str)]
    EXTRA_FIELDS = None

    def __init__(self, process_label, **kwargs):
        self.process_label = process_label
        self.work_chains_prompt = ipw.HTML(
            "<b>Select computed workflow or start a new one:</b>&nbsp;"
        )
        self.work_chains_selector = ipw.Dropdown(
            options=[("New workflow...", self._NO_PROCESS)],
            layout=ipw.Layout(min_width="300px", flex="1 1 auto"),
        )
        ipw.dlink(
            (self.work_chains_selector, "value"),
            (self, "value"),
            transform=lambda pk: None if pk is self._NO_PROCESS else pk,
        )

        if self.EXTRA_FIELDS is not None:
            fmt_extra = "".join([f"{{wc.{field[0]}}}" for field in self.EXTRA_FIELDS])
            self.FMT_WORKCHAIN += "\t" + fmt_extra

        self.refresh_work_chains_button = ipw.Button(description="Refresh")
        self.refresh_work_chains_button.on_click(self.refresh_work_chains)

        super().__init__(
            children=[
                self.work_chains_prompt,
                self.work_chains_selector,
                self.refresh_work_chains_button,
            ],
            **kwargs,
        )

        self.refresh_work_chains()

    def find_work_chains(self):
        builder = CalculationQueryBuilder()
        filters = builder.get_filters(
            process_label=self.process_label,
        )
        query_set = builder.get_query_set(
            filters=filters,
            order_by={"ctime": "desc"},
        )
        projections = ["pk", "ctime", "state"]
        projected = builder.get_projected(
            query_set,
            projections=projections,
        )

        for result in projected[1:]:
            process_info = dict(zip(projections, result))

            if self.EXTRA_FIELDS is not None:
                pk = process_info["pk"]
                extra_info = self.parse_extra_info(pk)

                yield make_dataclass("WorkChain", self.BASE_FIELDS + self.EXTRA_FIELDS)(
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
                    (self.FMT_WORKCHAIN.format(wc=wc), wc.pk)
                    for wc in self.find_work_chains()
                ]

                self.work_chains_selector.value = original_value
        finally:
            self.set_trait("busy", False)  # reenable the widget

    @tl.observe("value")
    def _observe_value(self, change):
        if change["old"] == change["new"]:
            return

        new = self._NO_PROCESS if change["new"] is None else change["new"]

        if new not in {pk for _, pk in self.work_chains_selector.options}:
            self.refresh_work_chains()

        self.work_chains_selector.value = new


class QeAppWorkChainSelector(WorkChainSelector):
    EXTRA_FIELDS = [("formula", str), ("relax_info", str), ("properties_info", str)]

    def __init__(self, **kwargs):
        super().__init__(process_label="QeAppWorkChain", **kwargs)

    def parse_extra_info(self, uuid):
        workchain = orm.load_node(uuid)
        formula = workchain.inputs.structure.get_formula()

        properties = []
        if "pdos" in workchain.inputs:
            properties.append("pdos")
        if "bands" in workchain.inputs:
            properties.append("bands")

        if "relax" in workchain.inputs:
            builder_parameters = workchain.get_extra("builder_parameters", {})
            relax_type = builder_parameters.get("relax_type")

            if relax_type != "none":
                relax_info = "structure is relaxed"
            else:
                relax_info = "static SCF calculation on structure"
        else:
            relax_info = "structure is not relaxed"

        if not properties:
            properties_info = ""
        else:
            properties_info = f"properties on {', '.join(properties)}"

        return {
            "formula": formula,
            "relax_info": relax_info,
            "properties_info": properties_info,
        }
