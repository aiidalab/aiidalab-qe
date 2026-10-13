import ipywidgets as ipw

from aiidalab_qe.common.panel import ConfigurationSettingsPanel

from .model import GeneralConfigurationSettingsModel


class GeneralConfigurationSettingsPanel(
    ConfigurationSettingsPanel[GeneralConfigurationSettingsModel],
):
    def render(self):
        if self.rendered:
            return

        self.clean_workdir = ipw.Checkbox(
            description="Delete the work directory after the calculation",
            indent=False,
            layout=ipw.Layout(width="fit-content", margin="5px 2px"),
        )
        ipw.link(
            (self._model, "clean_workdir"),
            (self.clean_workdir, "value"),
        )

        self.on_unhandled_failure = ipw.Dropdown(
            description="Action on failure:",
            options=[
                ("Abort", "abort"),
                ("Pause", "pause"),
                ("Restart once", "restart_once"),
                ("Restart and pause", "restart_and_pause"),
            ],
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "on_unhandled_failure"),
            (self.on_unhandled_failure, "value"),
        )

        # Total change setting
        self.total_charge = ipw.BoundedFloatText(
            min=-3,
            max=3,
            step=0.01,
            description="Total charge:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "total_charge"),
            (self.total_charge, "value"),
        )

        # Van der Waals setting
        self.van_der_waals = ipw.Dropdown(
            description="Van der Waals correction:",
            style={"description_width": "150px"},
        )
        ipw.dlink(
            (self._model, "van_der_waals_options"),
            (self.van_der_waals, "options"),
        )
        ipw.link(
            (self._model, "van_der_waals"),
            (self.van_der_waals, "value"),
        )

        self.children = [
            ipw.HTML(value="<h4>Workflow configuration</h4>"),
            self.clean_workdir,
            self.on_unhandled_failure,
            ipw.HTML(value="<h4>Global properties</h4>"),
            self.total_charge,
            self.van_der_waals,
        ]

        self.rendered = True
