# from aiidalab_qe.bands.result import Result
from aiidalab_qe.app.common.panel import OutlinePanel

from .result import Result
from .setting import Setting

# from .workchain import workchain_and_builder


class BandsOutline(OutlinePanel):
    title = "Electronic band structure"
    help = """The band structure workflow will
automatically detect the default path in reciprocal space using the
<a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">
SeeK-path tool</a>.
"""


bands = {
    "outline": BandsOutline,
    "setting": Setting,
    "result": Result,
    # "workchain": workchain_and_builder,
}
