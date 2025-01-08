import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.common.time import format_time, relative_time
from aiidalab_qe.common.widgets import TableWidget
from aiidalab_widgets_base.viewers import StructureDataViewer

from .model import StructureResultsModel


class StructureResultsPanel(ResultsPanel[StructureResultsModel]):
    def _render(self):
        if not hasattr(self, "widget"):
            structure = self._model.get_structure()
            self.widget = StructureDataViewer(structure=structure)
            # Select the Cell tab by default
            self.widget.configuration_box.selected_index = 2
            self.table_description = ipw.HTML("""
                <h4 style='margin: 10px 0;'>
                    Structure table information: Atom coordinates in Å
                </h4>
                <p style='margin: 5px 0; color: #555;'>
                    You can click on a row to select an atom. Multiple atoms
                    can be selected by clicking on additional rows. To unselect
                    an atom, click on the selected row again.
                </p>
            """)
            self.atom_coordinates_table = TableWidget()
            self._generate_table(structure.get_ase())

            # Basic widgets
            children = [
                self.widget,
                self.table_description,
                self.atom_coordinates_table,
            ]

            # Add structure info if it is a relaxed structure
            if "relax" in self._model.properties:
                self._initialize_structure_info(structure)
                children.insert(0, self.structure_info)

            # Add the children to the container
            self.results_container.children = tuple(children)
            self.atom_coordinates_table.observe(self._change_selection, "selected_rows")
            # Listen for changes in self.widget.displayed_selection and update the table
            self.widget.observe(self._update_table_selection, "displayed_selection")

        # HACK to resize the NGL viewer in cases where it auto-rendered when its
        # container was not displayed, which leads to a null width. This hack restores
        # the original dimensions.
        ngl = self.widget._viewer
        ngl._set_size("100%", "300px")
        ngl.control.zoom(0.0)

    def _initialize_structure_info(self, structure):
        self.structure_info = ipw.HTML(
            f"""
            <div style='line-height: 1.4;'>
                <strong>PK:</strong> {structure.pk}<br>
                <strong>Label:</strong> {structure.label}<br>
                <strong>Description:</strong> {structure.description}<br>
                <strong>Number of atoms:</strong> {len(structure.sites)}<br>
                <strong>Creation time:</strong> {format_time(structure.ctime)} ({relative_time(structure.ctime)})<br>
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
            # Format position values to two decimal places
            formatted_position = [f"{coord:.2f}" for coord in position]
            data.append([index, symbol, tag, *formatted_position])
        self.atom_coordinates_table.data = data

    def _change_selection(self, _):
        selected_indices = self.atom_coordinates_table.selected_rows
        self.widget.displayed_selection = selected_indices

    def _update_table_selection(self, change):
        selected_indices = change.new
        self.atom_coordinates_table.selected_rows = selected_indices
