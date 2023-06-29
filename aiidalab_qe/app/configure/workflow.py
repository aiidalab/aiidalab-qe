import ipywidgets as ipw

from aiidalab_qe.app.panel import Panel
from aiidalab_qe.app.utils import get_entry_items


class WorkChainSettings(Panel):
    identifier = "bands"

    relax_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Structure optimization</h4></div>"""
    )
    relax_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
        You have three options:<br>
        (1) Structure as is: perform a self consistent calculation using the structure provided as input.<br>
        (2) Atomic positions: perform a full relaxation of the internal atomic coordinates. <br>
        (3) Full geometry: perform a full relaxation for both the internal atomic coordinates and the cell vectors. </div>"""
    )
    properties_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Properties</h4></div>"""
    )

    def __init__(self, **kwargs):
        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "none"),
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ],
            value="positions_cell",
        )
        children = [
            self.relax_title,
            self.relax_help,
            self.relax_type,
            self.properties_title,
            ipw.HTML("Select which properties to calculate:"),
        ]
        self.properties = {}
        entries = get_entry_items("aiidalab_qe.property", "outline")
        for name, entry_point in entries.items():
            self.properties[name] = entry_point()
            children.append(self.properties[name])
        self.children = children
        super().__init__(
            **kwargs,
        )

    def get_panel_value(self):
        """Return the value of all the widgets in the panel as a dictionary.

        :return: a dictionary of the values of all the widgets in the panel.
        """
        parameters = {"relax_type": self.relax_type.value, "properties": {}}
        for name, property in self.properties.items():
            parameters["properties"][name] = property.run.value
        return parameters

    def set_panel_value(self, parameters):
        """Load a dictionary to set the value of the widgets in the panel.

        :param parameters: a dictionary of the values of all the widgets in the panel.
        """
        self.relax_type.value = parameters.get("relax_type", "positions_cell")
        for key, value in parameters.get("properties", {}).items():
            if key in self.properties:
                self.properties[key].run.value = value
