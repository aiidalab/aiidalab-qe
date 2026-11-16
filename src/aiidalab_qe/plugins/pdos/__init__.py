from aiidalab_qe.plugins import QeAppPlugin


class PdosPlugin(QeAppPlugin):
    @property
    def outline(self) -> str:
        return "Electronic projected density of states (PDOS)"

    def configuration(self):
        from .model import PdosConfigurationSettingsModel
        from .setting import PdosConfigurationSettingPanel

        return {
            "panel": PdosConfigurationSettingPanel,
            "model": PdosConfigurationSettingsModel,
        }

    def resources(self):
        from .resources import PdosResourceSettingsModel, PdosResourceSettingsPanel

        return {
            "panel": PdosResourceSettingsPanel,
            "model": PdosResourceSettingsModel,
        }

    def workchain(self):
        from .workchain import workchain_and_builder

        return workchain_and_builder

    def metadata(self):
        return {
            "process_labels": {
                "PdosWorkChain": "Projected density of states workflow",
                "DosCalculation": "Density of states calculation",
                "ProjwfcCalculation": "Orbital projection calculation",
            }
        }
