import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.common.widgets import TableWidget
from aiidalab_widgets_base.viewers import StructureDataViewer

from .model import StructureResultsModel


class StructureResultsPanel(ResultsPanel[StructureResultsModel]):
    def _render(self):
        if not hasattr(self, "widget"):
            structure = self._model.get_structure()
            self.widget = StructureDataViewer(structure=structure)
            self.widget.configuration_box.selected_index = (
                2  # Select the Cel tab by default
            )
            self.table_description = ipw.HTML(
                value="""
                <h4 style='margin: 10px 0;'>Structure Table Information: Atom Coordinates in Ångströms</h4>
                <p style='margin: 5px 0; color: #555;'>
                    You can click on a row to select an atom. Multiple atoms can be selected by clicking on additional rows. To unselect an atom, click on the selected row again.
                </p>
                """
            )
            self.atom_coordinates_table = TableWidget()
            self._generate_table(structure.get_ase())
            self.results_container.children = [
                self.widget,
                self.table_description,
                self.atom_coordinates_table,
            ]

        # HACK to resize the NGL viewer in cases where it auto-rendered when its
        # container was not displayed, which leads to a null width. This hack restores
        # the original dimensions.
        ngl = self.widget._viewer
        ngl._set_size("100%", "300px")
        ngl.control.zoom(0.0)

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
            # Format position values to two decimal places
            formatted_position = [f"{coord:.2f}" for coord in position]
            data.append([index, symbol, tag, *formatted_position])
        self.atom_coordinates_table.data = data
