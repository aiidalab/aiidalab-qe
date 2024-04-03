import ipywidgets as ipw
from aiidalab_qe.common.panel import ResultPanel


class Result(ResultPanel):
    title = "Bader Charge"
    workchain_labels = ["bader"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)
        self.summary_view = ipw.HTML()

    def _update_view(self):
        structure = self.node.inputs.bader.structure
        bader_charge = self.outputs.bader.bader.bader_charge.get_array("charge")
        self._generate_table(structure, bader_charge)
        self.children = [
            ipw.HBox(
                children=[self.summary_view],
                layout=ipw.Layout(justify_content="space-between", margin="10px"),
            ),
        ]

    def _generate_table(self, structure, bader_charge):
        # get index and element form AiiDA StructureData
        site_index = [site.kind_name for site in structure.sites]

        # Start of the HTML string for the table
        html_str = """<div class="custom-table" style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Bader Charge Table</h4>
                    <style>
                        .custom-table table, .custom-table th, .custom-table td {
                            border: 1px solid black;
                            border-collapse: collapse;
                            text-align: left;
                            padding: 8px;
                        }
                        .custom-table th, .custom-table td {
                            min-width: 120px;
                            word-wrap: break-word;
                        }
                        .custom-table table {
                            width: 100%;
                            font-size: 0.8em;
                        }
                    </style>
                    <table>
                        <tr>
                            <th>Site Index</th>
                            <th>Element</th>
                            <th>Bader Charge</th>
                        </tr>"""

        # Add rows to the table based on the bader_charge
        for i in range(len(site_index)):
            html_str += f"""
                        <tr>
                            <td>{i}</td>
                            <td>{site_index[i]}</td>
                            <td>{bader_charge[i]:1.3f}</td>
                        </tr>"""

        # Close the table and div tags
        html_str += """
                    </table>
                    </div>"""
        self.summary_view = ipw.HTML(html_str)
