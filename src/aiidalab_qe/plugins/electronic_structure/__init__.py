from aiidalab_qe.plugins import QeAppPlugin


class ElectronicStructurePlugin(QeAppPlugin):
    @property
    def outline(self) -> str:
        return "Electronic structure"

    def result(self):
        from .result import (
            ElectronicStructureResultsModel,
            ElectronicStructureResultsPanel,
        )

        return {
            "panel": ElectronicStructureResultsPanel,
            "model": ElectronicStructureResultsModel,
        }
