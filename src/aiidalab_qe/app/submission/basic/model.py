from __future__ import annotations

import os

import traitlets as tl

from aiida import orm
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.code import CodeModel, PwCodeModel
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import SettingsModel
from aiidalab_qe.common.widgets import (
    QEAppComputationalResourcesWidget,
)

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class BasicCodeModel(
    SettingsModel,
    HasInputStructure,
):
    """Model for the basic code setting."""

    dependencies = [
        "input_parameters",
        "input_structure",
    ]

    input_parameters = tl.Dict()

    codes = tl.Dict(
        key_trait=tl.Unicode(),  # code name
        value_trait=tl.Instance(CodeModel),  # code metadata
    )
    # this is a copy of the codes trait, which is used to trigger the update of the plugin
    basic_codes = tl.Dict(
        key_trait=tl.Unicode(),  # code name
        value_trait=tl.Instance(CodeModel),  # code metadata
    )

    plugin_mapping = tl.Dict(
        key_trait=tl.Unicode(),  # plugin identifier
        value_trait=tl.List(tl.Unicode()),  # list of code names
    )

    submission_blockers = tl.List(tl.Unicode())
    submission_warning_messages = tl.Unicode("")

    include = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Used by the code-setup thread to fetch code options
        # This is necessary to avoid passing the User object
        # between session in separate threads.
        self._default_user_email = orm.User.collection.get_default().email
        self._RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD = 10
        self._RUN_ON_LOCALHOST_VOLUME_WARN_THRESHOLD = 1000  # \AA^3

        self._ALERT_MESSAGE = """
            <div class="alert alert-{alert_class} alert-dismissible">
                <a href="#" class="close" data-dismiss="alert" aria-label="close">
                    &times;
                </a>
                <strong>{message}</strong>
            </div>
        """

    def refresh_codes(self):
        for _, code_model in self.codes.items():
            code_model.update(self._default_user_email)  # type: ignore

    def update_active_codes(self):
        for name, code_model in self.codes.items():
            if name != "quantumespresso.pw":
                code_model.deactivate()
        properties = self._get_properties()
        for identifier, code_names in self.plugin_mapping.items():
            if identifier in properties:
                for code_name in code_names:
                    self.codes[code_name].activate()

    def get_model_state(self):
        return {"codes": self.codes}

    def set_model_state(self, code_data: dict):
        for name, code_model in self.codes.items():
            if name in code_data:
                code_model.set_model_state(code_data[name])

    def add_code(self, identifier: str, code: CodeModel) -> CodeModel | None:
        """Add a code to the codes trait."""
        default_calc_job_plugin = code.default_calc_job_plugin
        name = default_calc_job_plugin.split(".")[-1]
        if default_calc_job_plugin not in self.codes:
            if default_calc_job_plugin == "quantumespresso.pw":
                code_model = PwCodeModel(
                    name=name,
                    description=name,
                    default_calc_job_plugin=default_calc_job_plugin,
                )
            else:
                code_model = CodeModel(
                    name=name,
                    description=name,
                    default_calc_job_plugin=default_calc_job_plugin,
                )
            self.codes[default_calc_job_plugin] = code_model
            return code_model
        # update the plugin mapping to keep track of which codes are associated with which plugins
        if identifier not in self.plugin_mapping:
            self.plugin_mapping[identifier] = [default_calc_job_plugin]
        else:
            self.plugin_mapping[identifier].append(default_calc_job_plugin)

    def get_code(self, name) -> CodeModel | None:
        if name in self.codes:  # type: ignore
            return self.codes[name]  # type: ignore

    def get_selected_codes(self) -> dict[str, dict]:
        return {
            name: code_model.get_model_state()
            for name, code_model in self.codes.items()
            if code_model.is_ready
        }

    def set_selected_codes(self, code_data=DEFAULT["codes"]):
        with self.hold_trait_notifications():
            for name, code_model in self.codes.items():
                if name in code_data:
                    code_model.set_model_state(code_data[name])

    def reset(self):
        """Reset the model to its default state."""
        for code_model in self.codes.values():
            code_model.reset()

    def _get_properties(self) -> list[str]:
        return self.input_parameters.get("workchain", {}).get("properties", [])

    def update_submission_blockers(self):
        self.submission_blockers = list(self._check_submission_blockers())

    def _check_submission_blockers(self):
        # No pw code selected (this is ignored while the setup process is running).
        pw_code = self._model.get_code("quantumespresso.pw")
        if pw_code and not pw_code.selected and not self.installing_qe:
            yield ("No pw code selected")

        # code related to the selected property is not installed
        properties = self._get_properties()
        message = "Calculating the {property} property requires code {code} to be set."
        for identifier, codes in self.get_model("basic").codes.items():
            if identifier in properties:
                for code in codes.values():
                    if not code.is_ready:
                        yield message.format(property=identifier, code=code.description)

        # check if the QEAppComputationalResourcesWidget is used
        for name, code in self.get_model("basic").codes.items():
            # skip if the code is not displayed, convenient for the plugin developer
            if not code.is_ready:
                continue
            if not issubclass(
                code.code_widget_class, QEAppComputationalResourcesWidget
            ):
                yield (
                    f"Error: hi, plugin developer, please use the QEAppComputationalResourcesWidget from aiidalab_qe.common.widgets for code {name}."
                )

    def check_resources(self):
        pw_code = self.get_code("quantumespresso.pw")

        if not self.input_structure or not pw_code.selected:
            return  # No code selected or no structure, so nothing to do

        num_cpus = pw_code.num_cpus * pw_code.num_nodes
        on_localhost = orm.load_node(pw_code.selected).computer.hostname == "localhost"
        num_sites = len(self.input_structure.sites)
        volume = self.input_structure.get_cell_volume()

        try:
            localhost_cpus = len(os.sched_getaffinity(0))
        except Exception:
            # Fallback, in some OS os.sched_getaffinity(0) is not supported
            # However, not so reliable in containers
            localhost_cpus = os.cpu_count()

        large_system = (
            num_sites > self._RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD
            or volume > self._RUN_ON_LOCALHOST_VOLUME_WARN_THRESHOLD
        )

        # Estimated number of CPUs for a run less than 12 hours.
        estimated_CPUs = self._estimate_min_cpus(num_sites, volume)

        # List of possible suggestions for warnings:
        suggestions = {
            "more_resources": f"<li>Increase the resources (total number of CPUs should be equal or more than {min(100,estimated_CPUs)}, if possible) </li>",
            "change_configuration": "<li>Review the configuration (e.g. choosing <i>fast protocol</i> - this will affect precision) </li>",
            "go_remote": "<li>Select a code that runs on a larger machine</li>",
            "avoid_overloading": "<li>Reduce the number of CPUs to avoid the overloading of the local machine </li>",
        }

        alert_message = ""
        if large_system and estimated_CPUs > num_cpus:
            # This part is in common between Warnings 1 (2):
            # (not) on localhost, big system and few cpus
            warnings_1_2 = (
                f"<span>&#9888;</span> Warning: The selected structure is large, with {num_sites} atoms "
                f"and a volume of {int(volume)} Å<sup>3</sup>, "
                "making it computationally demanding "
                "to run at the localhost. Consider the following: "
                if on_localhost
                else "to run in a reasonable amount of time. Consider the following: "
            )
            # Warning 1: on localhost, big system and few cpus
            alert_message += (
                f"{warnings_1_2}<ul>"
                + suggestions["more_resources"]
                + suggestions["change_configuration"]
                + "</ul>"
                if on_localhost
                else f"{warnings_1_2}<ul>"
                + suggestions["go_remote"]
                + suggestions["more_resources"]
                + suggestions["change_configuration"]
                + "</ul>"
            )
        if on_localhost and num_cpus / localhost_cpus > 0.8:
            # Warning-3: on localhost, more than half of the available cpus
            alert_message += (
                "<span>&#9888;</span> Warning: the selected pw.x code will run locally, but "
                f"the number of requested CPUs ({num_cpus}) is larger than the 80% of the available resources ({localhost_cpus}). "
                "Please be sure that your local "
                "environment has enough free CPUs for the calculation. Consider the following: "
                "<ul>"
                + suggestions["avoid_overloading"]
                + suggestions["go_remote"]
                + "</ul>"
            )
        print("alert_message: ", alert_message)

        self.submission_warning_messages = (
            ""
            if (on_localhost and num_cpus / localhost_cpus) <= 0.8
            and (not large_system or estimated_CPUs <= num_cpus)
            else self._ALERT_MESSAGE.format(
                alert_class="warning",
                message=alert_message,
            )
        )

    def _estimate_min_cpus(
        self,
        n,
        v,
        n0=9,
        v0=117,
        num_cpus0=4,
        t0=129.6,
        tmax=12 * 60 * 60,
        scf_cycles=5,
    ):
        """Estimate the minimum number of CPUs required to
        complete a task within a given time limit.

        Parameters
        ----------
        `n` : `int`
            The number of atoms in the system.
        `v` : `float`
            The volume of the system.
        `n0` : `int`, optional
            Reference number of atoms. Default is 9.
        `v0` : `float`, optional
            Reference volume. Default is 117.
        `num_cpus0` : `int`, optional
            Reference number of CPUs. Default is 4.
        `t0` : `float`, optional
            Reference time. Default is 129.6.
        `tmax` : `float`, optional
            Maximum time limit. Default is 12 hours.
        `scf_cycles` : `int`, optional
            Reference number of SCF cycles in a relaxation. Default is 5.

        Returns
        -------
        `int`
            The estimated minimum number of CPUs required.
        """
        import numpy as np

        return int(
            np.ceil(
                scf_cycles * num_cpus0 * (n / n0) ** 3 * (v / v0) ** 1.5 * t0 / tmax
            )
        )
