# from aiidalab_qe.bands.result import Result
from aiidalab_qe.app.submission.code.model import CodeModel
from aiidalab_qe.common.panel import PanelOutline

from .model import BandsModel
from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class BandsOutline(PanelOutline):
    title = "Electronic band structure"
    help = """The band structure workflow will
automatically detect the default path in reciprocal space using the
<a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">
SeeK-path tool</a>.
"""


bands = {
    "outline": BandsOutline,
    "model": BandsModel,
    "code": {
        "projwfc_bands": CodeModel(
            description="projwfc.x",
            default_calc_job_plugin="quantumespresso.projwfc",
        ),
    },
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
