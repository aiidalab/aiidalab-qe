"""Electronic structure results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.common.widgets import LoadingWidget

from .model import ElectronicStructureResultsModel


class ElectronicStructureResultsPanel(ResultsPanel[ElectronicStructureResultsModel]):
    def _render(self):
        apply_button = ipw.Button(
            description="Apply selection",
            button_style="primary",
            icon="pencil",
        )
        apply_button.on_click(self._render_property_results)

        property_selector = ipw.HBox(
            children=[apply_button],
            layout=ipw.Layout(grid_gap="10px"),
        )

        self.checkboxes: dict[str, ipw.Checkbox] = {}
        for identifier in self._model.identifiers:
            checkbox = ipw.Checkbox(
                description=self._model._TITLE_MAPPING[identifier],
                indent=False,
                value=self._model.has_partial_results(identifier),
                layout=ipw.Layout(width="fit-content"),
            )
            ipw.dlink(
                (self._model, "monitor_counter"),
                (checkbox, "disabled"),
                lambda _, cid=identifier: not self._model.has_partial_results(cid),
            )
            ipw.dlink(
                (checkbox, "value"),
                (apply_button, "disabled"),
                lambda _: not any(cb.value for cb in self.checkboxes.values()),
            )
            self.checkboxes[identifier] = checkbox
            property_selector.children += (checkbox,)

        self.sub_results_container = ipw.VBox()

        self.results_container.children = [
            ipw.HTML("Select one or more properties to plot:"),
            property_selector,
            self.sub_results_container,
        ]

    def _post_render(self):
        self._render_property_results()

    def _render_property_results(self, _=None):
        node_identifiers = [
            identifier
            for identifier, checkbox in self.checkboxes.items()
            if checkbox.value
        ]
        self._render_bands_pdos_widget(node_identifiers)

    def _render_bands_pdos_widget(self, node_identifiers):
        message = f"Loading {' + '.join(node_identifiers)} results"
        self.sub_results_container.children = [LoadingWidget(message)]
        nodes = {
            identifier: self._model.fetch_child_process_node(identifier)
            for identifier in node_identifiers
        }
        model = BandsPdosModel.from_nodes(**nodes)
        widget = BandsPdosWidget(model=model)
        widget.render()
        self.sub_results_container.children = [widget]
