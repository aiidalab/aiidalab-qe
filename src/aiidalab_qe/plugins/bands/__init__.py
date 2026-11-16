from aiidalab_qe.plugins import QeAppPlugin


class BandsPlugin(QeAppPlugin):
    @property
    def outline(self) -> str:
        return "Electronic band structure"

    def configuration(self):
        from .model import BandsConfigurationSettingsModel
        from .setting import BandsConfigurationSettingsPanel

        return {
            "panel": BandsConfigurationSettingsPanel,
            "model": BandsConfigurationSettingsModel,
        }

    def resources(self):
        from .resources import BandsResourceSettingsModel, BandsResourceSettingsPanel

        return {
            "panel": BandsResourceSettingsPanel,
            "model": BandsResourceSettingsModel,
        }

    def workchain(self):
        from .workchain import workchain_and_builder

        return workchain_and_builder

    def metadata(self):
        return {
            "process_labels": {
                "BandsWorkChain": "Electronic band structure workflow",
                "PwBandsWorkChain": "Electronic band structure workflow",
                "ProjwfcBandsWorkChain": "Electronic band structure workflow",
                "ProjwfcCalculation": "Orbital projection calculation",
            }
        }
