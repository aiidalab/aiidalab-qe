from aiidalab_qe.common.panel import OutlinePanel
from aiidalab_widgets_base import ComputationalResourcesWidget

from .result import Result
from .workchain import workchain_and_builder


class HpOutline(OutlinePanel):
    title = "HP calculation"
    help = """"""


hp_code = ComputationalResourcesWidget(
    description="hp.x",
    default_calc_job_plugin="quantumespresso.hp",
)


hp = {
    "outline": HpOutline,
    "code": {"hp": hp_code},
    "result": Result,
    "workchain": workchain_and_builder,
}
