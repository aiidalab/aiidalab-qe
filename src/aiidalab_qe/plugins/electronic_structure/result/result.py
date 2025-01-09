"""Electronic structure results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.common.widgets import LoadingWidget

from .model import ElectronicStructureResultsModel


class ElectronicStructureResultsPanel(ResultsPanel[ElectronicStructureResultsModel]):
    has_property_selector = False

    def _render(self):
        self.bands_pdos_container = ipw.VBox()
        children = []
        if self._model.needs_property_selector:
            children.append(self._render_property_selector())
            self.has_property_selector = True
        children.append(self.bands_pdos_container)
        self.results_container.children = children

    def _post_render(self):
        self._render_property_results()

    def _render_property_selector(self):
        apply_button = ipw.Button(
            description="Apply selection",
            button_style="primary",
            icon="pencil",
        )
        apply_button.on_click(self._render_property_results)

        property_selector = ipw.HBox(layout=ipw.Layout(grid_gap="10px"))

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

        property_selector.children += (apply_button,)

        return ipw.VBox(
            children=[
                ipw.HTML("""
                    <div>
                        <b>Select one or more properties to plot:</b>
                        <p>
                            You can choose to plot only bands, only PDOS, or both
                            combined in one plot. After making your selection, click
                            <span style="color: #2196f3;">
                                <i class="fa fa-pencil"></i> <b>Apply selection</b>
                            </span>
                            to proceed.
                        </p>
                    </div>
                """),
                property_selector,
            ]
        )

    def _render_property_results(self, _=None):
        node_identifiers = (
            [
                identifier
                for identifier, checkbox in self.checkboxes.items()
                if checkbox.value
            ]
            if self.has_property_selector
            else self._model.identifiers
        )
        self._render_bands_pdos_widget(node_identifiers)

    def _render_bands_pdos_widget(self, node_identifiers):
        message = f"Loading {' + '.join(node_identifiers)} results"
        self.bands_pdos_container.children = [LoadingWidget(message)]
        nodes = {
            **{
                identifier: self._model.fetch_child_process_node(identifier)
                for identifier in node_identifiers
            },
            "root": self._model.fetch_process_node(),
        }
        model = BandsPdosModel.from_nodes(**nodes)
        widget = BandsPdosWidget(model=model)
        widget.render()
        self.bands_pdos_container.children = [widget]
