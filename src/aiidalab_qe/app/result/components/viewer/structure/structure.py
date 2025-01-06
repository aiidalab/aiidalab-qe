import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel
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
            self.atom_coordinates_table = self.create_structure_table(
                structure.get_ase()
            )
            self.results_container.children = [self.widget, self.atom_coordinates_table]

        # HACK to resize the NGL viewer in cases where it auto-rendered when its
        # container was not displayed, which leads to a null width. This hack restores
        # the original dimensions.
        ngl = self.widget._viewer
        ngl._set_size("100%", "300px")
        ngl.control.zoom(0.0)

    def create_structure_table(self, structure):
        # Extract data from ase structure
        positions = structure.positions
        chemical_symbols = structure.get_chemical_symbols()
        tags = structure.get_tags()

        # Define styles for table and cells
        table_style = "border: 1px solid black; border-collapse: collapse; text-align: center; width: auto; margin: 10px;"
        cell_style = "border: 1px solid black; padding: 4px; text-align: center;"

        # Start table HTML with headers
        headers = ["Atom<br>index", "Chemical<br>symbol", "Tag", "x", "y", "z"]
        html_content = f"<table style='{table_style}'><tr>"
        html_content += "".join(
            [f"<th style='{cell_style}'>{header}</th>" for header in headers]
        )
        html_content += "</tr>"

        # Populate the table rows
        for index, (symbol, tag, position) in enumerate(
            zip(chemical_symbols, tags, positions), start=1
        ):
            html_content += f"""
            <tr>
                <td style="{cell_style}">{index}</td>
                <td style="{cell_style}">{symbol}</td>
                <td style="{cell_style}">{tag}</td>
                <td style="{cell_style}">{position[0]:.3f}</td>
                <td style="{cell_style}">{position[1]:.3f}</td>
                <td style="{cell_style}">{position[2]:.3f}</td>
            </tr>
            """

        # Finish the table
        html_content += "</table>"

        # Create an HTML widget
        html_widget = ipw.HTML(value=html_content)
        return html_widget
