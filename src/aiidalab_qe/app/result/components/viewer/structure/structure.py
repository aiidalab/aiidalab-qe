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

        html_content = """
        <table style="border: 1px solid black; border-collapse: collapse; text-align: center; width: auto; margin: 10px;">
        <tr>
            <th style="border: 1px solid black; padding: 4px; text-align: center;">Atom<br>index</th>
            <th style="border: 1px solid black; padding: 4px; text-align: center;">Chemical<br>symbol</th>
            <th style="border: 1px solid black; padding: 4px; text-align: center;">Tag</th>
            <th style="border: 1px solid black; padding: 4px; text-align: center;">x</th>
            <th style="border: 1px solid black; padding: 4px; text-align: center;">y</th>
            <th style="border: 1px solid black; padding: 4px; text-align: center;">z</th>
        </tr>
        """

        # populate the table rows
        for index, (symbol, tag, position) in enumerate(
            zip(chemical_symbols, tags, positions)
        ):
            html_content += f"""
            <tr>
            <td style="border: 1px solid black; padding: 4px;">{index+1}</td>
            <td style="border: 1px solid black; padding: 4px;">{symbol}</td>
            <td style="border: 1px solid black; padding: 4px;">{tag}</td>
            <td style="border: 1px solid black; padding: 4px;">{position[0]:.3f}</td>
            <td style="border: 1px solid black; padding: 4px;">{position[1]:.3f}</td>
            <td style="border: 1px solid black; padding: 4px;">{position[2]:.3f}</td>
            </tr>
            """

        # finish the table
        html_content += "</table>"

        # Create an HTML widget
        html_widget = ipw.HTML(value=html_content)

        return html_widget
