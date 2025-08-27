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
            self.clean_workdir,
            self.total_charge,
            self.van_der_waals,
        ]

        self.rendered = True
