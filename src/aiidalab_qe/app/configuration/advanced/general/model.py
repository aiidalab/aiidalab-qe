from __future__ import annotations

import typing as t

import traitlets as tl

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS

from ..subsettings import AdvancedCalculationSubSettingsModel

DEFAULT = t.cast(dict, DEFAULT_PARAMETERS)


class GeneralConfigurationSettingsModel(AdvancedCalculationSubSettingsModel):
    title = "General"
    identifier = "general"

    clean_workdir = tl.Bool(DEFAULT["advanced"]["clean_workdir"])
    total_charge = tl.Float(DEFAULT["advanced"]["tot_charge"])
    van_der_waals_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["None", "none"],
            ["Grimme-D3", "dft-d3"],
            ["Grimme-D3BJ", "dft-d3bj"],
            ["Grimme-D3M", "dft-d3m"],
            ["Grimme-D3MBJ", "dft-d3mbj"],
            ["Tkatchenko-Scheffler", "ts-vdw"],
        ],
    )
    van_der_waals = tl.Unicode(DEFAULT["advanced"]["vdw_corr"])

    dftd3_version = {
        "dft-d3": 3,
        "dft-d3bj": 4,
        "dft-d3m": 5,
        "dft-d3mbj": 6,
    }

    def reset(self):
        with self.hold_trait_notifications():
            self.total_charge = self._get_default("total_charge")
            self.van_der_waals = self._get_default("van_der_waals")
