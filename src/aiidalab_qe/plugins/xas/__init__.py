from importlib import resources

import yaml
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget

from aiidalab_qe.common.panel import OutlinePanel
from aiidalab_qe.plugins import xas as xas_folder

from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder

PSEUDO_TOC = yaml.safe_load(resources.read_text(xas_folder, "pseudo_toc.yaml"))


class XasOutline(OutlinePanel):
    title = "X-ray absorption spectroscopy (XAS)"
    help = """"""


xs_code = QEAppComputationalResourcesWidget(
    description="xspectra.x", default_calc_job_plugin="quantumespresso.xspectra"
)

xas = {
    "outline": XasOutline,
    "code": {"xspectra": xs_code},
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
