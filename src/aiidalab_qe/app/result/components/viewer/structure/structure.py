import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.common.widgets import TableWidget
from aiidalab_widgets_base.viewers import StructureDataViewer

from .model import StructureResultsModel


class StructureResultsPanel(ResultsPanel[StructureResultsModel]):
    def __init__(self, model: StructureResultsModel, **kwargs):
        super().__init__(model, **kwargs)
        self._model.observe(
            self._on_selected_view_change,
            "selected_view",
        )

    def _render(self):
        if hasattr(self, "widget"):
            # HACK to resize the NGL viewer in cases where it auto-rendered when its
            # container was not displayed, which leads to a null width. This hack
            # restores the original dimensions.
            ngl = self.widget._viewer
            ngl._set_size("100%", "300px")
            ngl.control.zoom(0.0)
            return

        self.header = ipw.HTML()
        ipw.dlink(
            (self._model, "header"),
            (self.header, "value"),
        )

        self.sub_header = ipw.HTML()
        ipw.dlink(
            (self._model, "sub_header"),
            (self.sub_header, "value"),
        )

        self.view_toggle_button = ipw.Button(
            icon="eye",
            layout=ipw.Layout(
                display="block" if self._model.is_relaxed else "none",
                width="125px",
            ),
        )
        ipw.dlink(
            (self._model, "selected_view"),
            (self.view_toggle_button, "description"),
            lambda view: f"View {'initial' if view == 'relaxed' else 'relaxed'}",
        )
        ipw.dlink(
            (self._model, "monitor_counter"),
            (self.view_toggle_button, "disabled"),
            lambda _: not self._model.has_results,
        )
        ipw.dlink(
            (self.view_toggle_button, "disabled"),
            (self.view_toggle_button, "tooltip"),
            lambda disabled: "Waiting for results"
            if disabled
            else "Toggle between the initial and relaxed structures",
        )

        self.view_toggle_button.on_click(self._toggle_view)

        self.structure_info = ipw.HTML(layout=ipw.Layout(margin="0"))
        ipw.dlink(
            (self._model, "info"),
            (self.structure_info, "value"),
        )

        self.header_box = ipw.HBox(
            children=[
                ipw.VBox(
                    children=[
                        ipw.VBox(
                            children=[
                                self.header,
                                self.sub_header,
                            ],
                        ),
                        self.view_toggle_button,
                    ],
                    layout=ipw.Layout(justify_content="space-between"),
                ),
                ipw.VBox(
                    children=[
                        self.structure_info,
                    ],
                    layout=ipw.Layout(justify_content="flex-end"),
                ),
            ],
            layout=ipw.Layout(grid_gap="1em"),
        )

        self.widget = StructureDataViewer()
        ipw.dlink(
            (self._model, "structure"),
            (self.widget, "structure"),
        )

        self.widget.configuration_box.selected_index = 2  # select the Cell tab

        self.atom_coordinates_table = TableWidget()
        self.atom_coordinates_table.add_class("atom-coordinates-table")
        ipw.dlink(
            (self._model, "table_data"),
            (self.atom_coordinates_table, "data"),
        )

        ipw.link(
            (self.widget, "displayed_selection"),
            (self.atom_coordinates_table, "selected_rows"),
        )

        self.results_container.children = [
            self.header_box,
            self.widget,
            ipw.HTML("""
                <h4 style='margin: 10px 0;'>
                    Structure information: Atom coordinates in Å
                </h4>
                <p style='margin: 5px 0; color: #555;'>
                    You can click on a row to select an atom. Multiple atoms can be
                    selected by clicking on additional rows. To unselect an atom, click
                    on the selected row again.
                </p>
            """),
            self.atom_coordinates_table,
        ]
        # same reason as above for the hack
        self.widget._viewer.handle_resize()
        # trigger the view update
        structure = self._model.structure
        self._model.structure = None
        self._model.structure = structure

    def _on_process_change(self, _):
        super()._on_process_change(_)
        if self.rendered:
            self.view_toggle_button.layout.display = (
                "block" if self._model.is_relaxed else "none"
            )

    def _on_selected_view_change(self, _):
        self._model.update()

    def _toggle_view(self, _):
        self._model.toggle_selected_view()
