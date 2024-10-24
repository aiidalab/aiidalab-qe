import traitlets as tl

from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget


class CodeModel(tl.HasTraits):
    is_active = tl.Bool(False)
    selected = tl.Bool(False)
    options = tl.List(
        trait=tl.Tuple(tl.Unicode(), tl.Unicode()),  # code option (label, uuid)
        default_value=[],
    )
    parameters = tl.Dict()

    def __init__(
        self,
        *,
        name="pw",
        description,
        default_calc_job_plugin,
        setup_widget_class=QEAppComputationalResourcesWidget,
    ):
        self.name = name
        self.description = description
        self.default_calc_job_plugin = default_calc_job_plugin
        self.setup_widget_class = setup_widget_class
        self.is_rendered = False
        self.is_loaded = False

    @property
    def is_ready(self):
        return self.is_loaded and self.is_active and self.selected

    def activate(self):
        self.is_active = True

    def deactivate(self):
        self.is_active = False

    def get_setup_widget(self) -> QEAppComputationalResourcesWidget:
        if not self.is_loaded:
            self._setup_widget = self.setup_widget_class(
                description=self.description,
                default_calc_job_plugin=self.default_calc_job_plugin,
            )
            self.is_loaded = True
        return self._setup_widget


CodesDict = dict[str, CodeModel]
PluginCodes = dict[str, CodesDict]
