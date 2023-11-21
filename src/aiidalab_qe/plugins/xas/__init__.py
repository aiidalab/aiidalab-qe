from aiidalab_widgets_base import ComputationalResourcesWidget

from aiidalab_qe.common.panel import OutlinePanel

from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class XasOutline(OutlinePanel):
    title = "X-ray absorption spectroscopy (XAS)"
    help = """"""


xs_code = ComputationalResourcesWidget(
    description="xspectra.x", default_calc_job_plugin="quantumespresso.xspectra"
)

xas = {
    "outline": XasOutline,
    "code": {"xspectra": xs_code},
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
