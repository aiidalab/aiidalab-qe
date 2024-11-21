import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.common import NotExistent
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import (
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)


class CodeModel(Model):
    is_active = tl.Bool(False)
    options = tl.List(
        trait=tl.Tuple(tl.Unicode(), tl.Unicode()),  # code option (label, uuid)
        default_value=[],
    )
    selected = tl.Unicode(default_value=None, allow_none=True)
    num_nodes = tl.Int(1)
    num_cpus = tl.Int(1)
    ntasks_per_node = tl.Int(1)
    cpus_per_task = tl.Int(1)
    max_wallclock_seconds = tl.Int(3600 * 12)
    allow_hidden_codes = tl.Bool(False)
    allow_disabled_computers = tl.Bool(False)

    def __init__(
        self,
        *,
        name="",
        description,
        default_calc_job_plugin,
        code_widget_class=QEAppComputationalResourcesWidget,
    ):
        self.name = name
        self.description = description
        self.default_calc_job_plugin = default_calc_job_plugin
        self.code_widget_class = code_widget_class
        self.is_rendered = False

        ipw.dlink(
            (self, "num_cpus"),
            (self, "ntasks_per_node"),
        )

    @property
    def is_ready(self):
        return self.is_active and bool(self.selected)

    def activate(self):
        self.is_active = True

    def deactivate(self):
        self.is_active = False

    def update(self, user_email: str):
        if not self.options:
            self.options = self._get_codes(user_email)
            self.selected = self.options[0][1] if self.options else None

    def get_model_state(self) -> dict:
        return {
            "code": self.selected,
            "nodes": self.num_nodes,
            "cpus": self.num_cpus,
            "ntasks_per_node": self.ntasks_per_node,
            "cpus_per_task": self.cpus_per_task,
            "max_wallclock_seconds": self.max_wallclock_seconds,
        }

    def set_model_state(self, parameters):
        self.selected = self._get_uuid(parameters["code"])
        self.num_nodes = parameters.get("nodes", 1)
        self.num_cpus = parameters.get("cpus", 1)
        self.ntasks_per_node = parameters.get("ntasks_per_node", 1)
        self.cpus_per_task = parameters.get("cpus_per_task", 1)
        self.max_wallclock_seconds = parameters.get("max_wallclock_seconds", 3600 * 12)

    def _get_uuid(self, identifier):
        if not self.selected:
            try:
                uuid = orm.load_code(identifier).uuid
            except NotExistent:
                uuid = None
            # If the code was imported from another user, it is not usable
            # in the app and thus will not be considered as an option!
            self.selected = uuid if uuid in [opt[1] for opt in self.options] else None
        return self.selected

    def _get_codes(self, user_email: str):
        # set default user_email if not provided
        user_email = user_email or orm.User.collection.get_default().email
        user = orm.User.collection.get(email=user_email)

        filters = (
            {"attributes.input_plugin": self.default_calc_job_plugin}
            if self.default_calc_job_plugin
            else {}
        )

        codes = (
            orm.QueryBuilder()
            .append(
                orm.Code,
                filters=filters,
            )
            .all(flat=True)
        )

        return [
            (self._full_code_label(code), code.uuid)
            for code in codes
            if code.computer.is_user_configured(user)
            and (self.allow_hidden_codes or not code.is_hidden)
            and (self.allow_disabled_computers or code.computer.is_user_enabled(user))
        ]

    @staticmethod
    def _full_code_label(code):
        return f"{code.label}@{code.computer.label}"


class PwCodeModel(CodeModel):
    override = tl.Bool(False)
    npool = tl.Int(1)

    def __init__(
        self,
        *,
        name="",
        description,
        default_calc_job_plugin,
        code_widget_class=PwCodeResourceSetupWidget,
    ):
        super().__init__(
            name=name,
            description=description,
            default_calc_job_plugin=default_calc_job_plugin,
            code_widget_class=code_widget_class,
        )

    def get_model_state(self) -> dict:
        parameters = super().get_model_state()
        parameters["parallelization"] = {"npool": self.npool} if self.override else {}
        return parameters

    def set_model_state(self, parameters):
        super().set_model_state(parameters)
        if "parallelization" in parameters and "npool" in parameters["parallelization"]:
            self.override = True
            self.npool = parameters["parallelization"].get("npool", 1)


CodesDict = dict[str, CodeModel]
PluginCodes = dict[str, CodesDict]
