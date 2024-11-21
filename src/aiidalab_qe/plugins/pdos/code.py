"""Panel for pdos plugin."""

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiidalab_qe.common.code.model import CodeModel, PwCodeModel
from aiidalab_qe.common.panel import SettingsModel, SettingsPanel
from aiidalab_qe.common.widgets import (
    LoadingWidget,
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)


class PdosCodeModel(SettingsModel):
    """Model for the pdos structure plugin."""

    dependencies = []

    codes = {
        "pw": PwCodeModel(
            name="pw.x",
            description="pw.x",
            default_calc_job_plugin="quantumespresso.pw",
        ),
        "dos": CodeModel(
            name="dos.x",
            description="dos.x",
            default_calc_job_plugin="quantumespresso.dos",
        ),
        "projwfc": CodeModel(
            name="projwfc.x",
            description="projwfc.x",
            default_calc_job_plugin="quantumespresso.projwfc",
        ),
    }
    global_codes = tl.Dict(
        key_trait=tl.Unicode(),  # code name
        value_trait=tl.Dict(),  # code metadata
    )

    submission_blockers = tl.List(tl.Unicode())
    submission_warning_messages = tl.Unicode("")
    override = tl.Bool(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Used by the code-setup thread to fetch code options
        # This is necessary to avoid passing the User object
        # between session in separate threads.
        self._default_user_email = orm.User.collection.get_default().email

    def refresh_codes(self):
        for _, code_model in self.codes.items():
            code_model.update(self._default_user_email)  # type: ignore

    def update_code_from_global(self):
        # skip the sync if the user has overridden the settings
        if self.override:
            return
        for _, code_model in self.codes.items():
            default_calc_job_plugin = code_model.default_calc_job_plugin
            if default_calc_job_plugin in self.global_codes:
                code_data = self.global_codes[default_calc_job_plugin]
                code_model.set_model_state(code_data)

    def get_model_state(self):
        codes = {name: model.get_model_state() for name, model in self.codes.items()}
        return {
            "codes": codes,
            "override": self.override,
        }

    def set_model_state(self, code_data: dict):
        for name, code_model in self.codes.items():
            if name in code_data:
                code_model.set_model_state(code_data[name])

    def reset(self):
        """Reset the model to its default state."""
        for code_model in self.codes.values():
            code_model.reset()


class PdosCodeSettings(SettingsPanel[PdosCodeModel]):
    title = "PDOS"
    identifier = "pdos"

    def render(self):
        if self.rendered:
            return
        self.code_widgets_container = ipw.VBox()
        self.code_widgets = {}
        self.children = [
            self.code_widgets_container,
        ]

        self.rendered = True

        for code_model in self._model.codes.values():
            self._toggle_code(code_model)
        return self.code_widgets_container

    def _toggle_code(self, code_model: CodeModel):
        if not self.rendered:
            return
        if not code_model.is_rendered:
            loading_message = LoadingWidget(f"Loading {code_model.name} code")
            self.code_widgets_container.children += (loading_message,)
        if code_model.name not in self.code_widgets:
            code_widget = code_model.code_widget_class(
                description=code_model.description,
                default_calc_job_plugin=code_model.default_calc_job_plugin,
            )
            self.code_widgets[code_model.name] = code_widget
        else:
            code_widget = self.code_widgets[code_model.name]
        # TODO: hide code if it is not used
        # code_widget.layout.display = "block" if code_model.is_active else "none"
        if not code_model.is_rendered:
            self._render_code_widget(code_model, code_widget)

    def _render_code_widget(
        self,
        code_model: CodeModel,
        code_widget: QEAppComputationalResourcesWidget,
    ):
        code_model.update()
        ipw.dlink(
            (code_model, "options"),
            (code_widget.code_selection.code_select_dropdown, "options"),
        )
        ipw.link(
            (code_model, "selected"),
            (code_widget.code_selection.code_select_dropdown, "value"),
        )
        ipw.dlink(
            (code_model, "selected"),
            (code_widget.code_selection.code_select_dropdown, "disabled"),
            lambda selected: not selected,
        )
        ipw.link(
            (code_model, "num_cpus"),
            (code_widget.num_cpus, "value"),
        )
        ipw.link(
            (code_model, "num_nodes"),
            (code_widget.num_nodes, "value"),
        )
        ipw.link(
            (code_model, "ntasks_per_node"),
            (code_widget.resource_detail.ntasks_per_node, "value"),
        )
        ipw.link(
            (code_model, "cpus_per_task"),
            (code_widget.resource_detail.cpus_per_task, "value"),
        )
        ipw.link(
            (code_model, "max_wallclock_seconds"),
            (code_widget.resource_detail.max_wallclock_seconds, "value"),
        )
        if isinstance(code_widget, PwCodeResourceSetupWidget):
            ipw.link(
                (code_model, "override"),
                (code_widget.parallelization.override, "value"),
            )
            ipw.link(
                (code_model, "npool"),
                (code_widget.parallelization.npool, "value"),
            )
            code_model.observe(
                self._on_pw_code_resource_change,
                ["num_cpus", "num_nodes"],
            )
        code_widgets = self.code_widgets_container.children[:-1]  # type: ignore

        self.code_widgets_container.children = [*code_widgets, code_widget]
        code_model.is_rendered = True
