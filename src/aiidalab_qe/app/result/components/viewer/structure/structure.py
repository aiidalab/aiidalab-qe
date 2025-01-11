import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.common.time import format_time, relative_time
from aiidalab_qe.common.widgets import TableWidget
from aiidalab_widgets_base.viewers import StructureDataViewer

from .model import StructureResultsModel


class StructureResultsPanel(ResultsPanel[StructureResultsModel]):
    def _render(self):
        if hasattr(self, "widget"):
            # HACK to resize the NGL viewer in cases where it auto-rendered when its
            # container was not displayed, which leads to a null width. This hack restores
            # the original dimensions.
            ngl = self.widget._viewer
            ngl._set_size("100%", "300px")
            ngl.control.zoom(0.0)
            return

        structure = self._model.get_structure()
        self.widget = StructureDataViewer(structure=structure)

        self.widget.configuration_box.selected_index = 2  # select the Cell tab

        self.atom_coordinates_table = TableWidget()
        self.atom_coordinates_table.add_class("atom-coordinates-table")

        self._generate_table(structure.get_ase())

        structure_info = self._get_structure_info(structure)

        ipw.link(
            (self.widget, "displayed_selection"),
            (self.atom_coordinates_table, "selected_rows"),
        )

        self.results_container.children = [
            structure_info,
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

    def _get_structure_info(self, structure):
        formatted = format_time(structure.ctime)
        relative = relative_time(structure.ctime)
        return ipw.HTML(
            f"""
            <div style='line-height: 1.4;'>
                <strong>PK:</strong> {structure.pk}<br>
                <strong>Label:</strong> {structure.label}<br>
                <strong>Description:</strong> {structure.description}<br>
                <strong>Number of atoms:</strong> {len(structure.sites)}<br>
                <strong>Creation time:</strong> {formatted} ({relative})<br>
            </div>
            """
        )

    def _generate_table(self, structure):
        data = [
            [
                "Atom index",
                "Chemical symbol",
                "Tag",
                "x (Å)",
                "y (Å)",
                "z (Å)",
            ]
        ]
        positions = structure.positions
        chemical_symbols = structure.get_chemical_symbols()
        tags = structure.get_tags()

        for index, (symbol, tag, position) in enumerate(
            zip(chemical_symbols, tags, positions), start=1
        ):
            formatted_position = [f"{coord:.2f}" for coord in position]
            data.append([index, symbol, tag, *formatted_position])

        self.atom_coordinates_table.data = data
