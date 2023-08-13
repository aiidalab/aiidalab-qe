"""Electronic structure results view widgets

"""
import ipywidgets as ipw

from aiidalab_qe.app.common.panel import ResultPanel


def export_data(work_chain_node, group_dos_by="atom"):
    from aiidalab_qe.app.plugins.bands.result import export_bands_data
    from aiidalab_qe.app.plugins.pdos.result import export_pdos_data

    if "pdos" in work_chain_node.outputs:
        dos = export_pdos_data(work_chain_node.outputs.pdos, group_dos_by=group_dos_by)
        fermi_energy = dos["fermi_energy"] if dos else None
    else:
        dos = None
        fermi_energy = None
    if "bands" in work_chain_node.outputs:
        bands = export_bands_data(work_chain_node.outputs.bands, fermi_energy)
    else:
        bands = None

    return dict(
        bands=bands,
        dos=dos,
    )


class ElectronicStructure(ResultPanel):
    title = "Electronic Structure"

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, identifier="electronic", **kwargs)

    def _update_view(self):
        from widget_bandsplot import BandsPlotWidget

        group_dos_by = ipw.ToggleButtons(
            options=[
                ("Atom", "atom"),
                ("Orbital", "angular"),
            ],
            value="atom",
        )
        settings = ipw.VBox(
            children=[
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "DOS grouped by:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="120px"
                            ),
                        ),
                        group_dos_by,
                    ]
                ),
            ],
            layout={"margin": "0 0 30px 30px"},
        )
        #
        data = export_data(self.node, group_dos_by=group_dos_by.value)
        bands_data = data.get("bands", None)
        dos_data = data.get("dos", None)
        if bands_data is None and dos_data is None:
            return
        _bands_plot_view = BandsPlotWidget(
            bands=bands_data,
            dos=dos_data,
        )

        def response(change):
            data = export_data(self.node, group_dos_by=group_dos_by.value)
            bands_data = data.get("bands", None)
            dos_data = data.get("dos", None)
            if bands_data is None and dos_data is None:
                return
            _bands_plot_view = BandsPlotWidget(
                bands=bands_data,
                dos=dos_data,
            )
            self.children = [
                settings,
                _bands_plot_view,
            ]

        group_dos_by.observe(response, names="value")
        # update the electronic structure tab
        self.children = [
            settings,
            _bands_plot_view,
        ]
