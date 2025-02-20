from __future__ import annotations

import traitlets as tl

from aiida import orm
from aiidalab_qe.common.code import CodeModel, PwCodeModel
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import ResourceSettingsModel
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget


class GlobalResourceSettingsModel(
    ResourceSettingsModel,
    HasInputStructure,
):
    """Model for the global code setting."""

    title = "Global resources"
    identifier = "global"

    dependencies = [
        "input_structure",
        "input_parameters",
    ]

    input_parameters = tl.Dict()

    plugin_overrides = tl.List(tl.Unicode())
    plugin_overrides_notification = tl.Unicode("")

    include = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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

        self.plugin_mapping: dict[str, list[str]] = {}

    def update_global_codes(self):
        self.global_codes = self.get_model_state()["codes"]

    def update_active_codes(self):
        for identifier, code_model in self.get_models():
            if identifier != "quantumespresso__pw":
                code_model.deactivate()
        properties = self._get_properties()
        for identifier, code_names in self.plugin_mapping.items():
            if identifier in properties:
                for code_name in code_names:
                    self.get_model(code_name).activate()

    def update_plugin_overrides_notification(self):
        if self.plugin_overrides:
            formatted = "\n".join(
                f"<li>{plugin}</li>" for plugin in self.plugin_overrides
            )
            self.plugin_overrides_notification = f"""
                <div class="alert alert-info" style="margin-top: 5px; margin-bottom: 0">
                    <strong>Currently overriding computational resources for:</strong>
                    <ul>
                        {formatted}
                    </ul>
                </div>
            """
        else:
            self.plugin_overrides_notification = ""

    def add_global_model(
        self,
        identifier: str,
        code_model: CodeModel,
    ) -> CodeModel | None:
        """Registers a code with this model."""
        base_code_model = None
        default_calc_job_plugin = code_model.default_calc_job_plugin
        name = default_calc_job_plugin.split(".")[-1]
        # "." in the model key means nested models
        model_key = default_calc_job_plugin.replace(".", "__")

        if not self.has_model(default_calc_job_plugin):
            if default_calc_job_plugin == "quantumespresso.pw":
                base_code_model = PwCodeModel(
                    name=name,
                    description=name,
                    default_calc_job_plugin=default_calc_job_plugin,
                )
                base_code_model.activate()
            else:
                base_code_model = CodeModel(
                    name=name,
                    description=name,
                    default_calc_job_plugin=default_calc_job_plugin,
                )
            self.add_model(model_key, base_code_model)

        if identifier not in self.plugin_mapping:
            self.plugin_mapping[identifier] = [model_key]
        else:
            self.plugin_mapping[identifier].append(model_key)

        return base_code_model

    def check_resources(self):
        if not self.has_model("quantumespresso__pw"):
            return

        pw_code_model = self.get_model("quantumespresso__pw")
        protocol = self.input_parameters.get("workchain", {}).get("protocol", "fast")

        if not self.input_structure or not pw_code_model.selected:
            return  # No code selected or no structure, so nothing to do

        num_cpus = pw_code_model.num_cpus * pw_code_model.num_nodes
        on_localhost = (
            orm.load_node(pw_code_model.selected).computer.hostname == "localhost"
        )
        num_sites = len(self.input_structure.sites)
        volume = self.input_structure.get_cell_volume()

        code = orm.load_node(pw_code_model.selected)
        machine_cpus = code.computer.get_default_mpiprocs_per_machine() or 1

        large_system = (
            num_sites > self._RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD
            or volume > self._RUN_ON_LOCALHOST_VOLUME_WARN_THRESHOLD
        )

        # Estimated number of CPUs for a run less than 12 hours.
        estimated_CPUs = self._estimate_min_cpus(num_sites, volume)

        # List of possible suggestions for warnings:
        suggestions = {
            "more_resources": f"<li>Increase the resources (total number of CPUs should be equal or more than {min(100, estimated_CPUs)}, if possible) </li>",
            "change_configuration": "<li>Review the configuration (e.g. choosing <i>fast protocol</i> - this will affect precision) </li>",
            "go_remote": "<li>Select a code that runs on a larger machine</li>",
            "avoid_overloading": "<li>Reduce the number of CPUs to avoid the overloading of the local machine </li>",
        }

        alert_message = ""
        if (
            on_localhost and protocol == "precise" and machine_cpus < 4
        ):  # This means that we are in a small deployment.
            alert_message += (
                f"Warning: The selected protocol is {protocol}, which is computationally demanding to run on localhost. "
                f"Consider the following: <ul>"
                + suggestions["change_configuration"]
                + "</ul>"
            )
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
        if on_localhost and num_cpus / machine_cpus > 0.8:
            # Warning-3: on localhost, more than half of the available cpus
            alert_message += (
                "<span>&#9888;</span> Warning: the selected pw.x code will run locally, but "
                f"the number of requested CPUs ({num_cpus}) is larger than the 80% of the available resources ({machine_cpus}). "
                "Please be sure that your local "
                "environment has enough free CPUs for the calculation. Consider the following: "
                "<ul>"
                + suggestions["avoid_overloading"]
                + suggestions["go_remote"]
                + "</ul>"
            )

        self.warning_messages = (
            ""
            if (on_localhost and num_cpus / machine_cpus) <= 0.8
            and (not large_system or estimated_CPUs <= num_cpus)
            else self._ALERT_MESSAGE.format(
                alert_class="warning",
                message=alert_message,
            )
        )

    def _get_properties(self) -> list[str]:
        return self.input_parameters.get("workchain", {}).get("properties", [])

    def _check_blockers(self):
        # No pw code selected
        pw_code_model = self.get_model("quantumespresso__pw")
        if pw_code_model and not pw_code_model.selected:
            yield ("No pw code selected")

        # Code related to the selected property is not installed
        properties = self._get_properties()
        message = "Calculating the {property} property requires code {code} to be set."
        for identifier, code_names in self.plugin_mapping.items():
            if identifier in properties:
                for name in code_names:
                    code_model = self.get_model(name)
                    if not code_model.is_ready:
                        yield message.format(
                            property=name,
                            code=code_model.description,
                        )

        # Check if the QEAppComputationalResourcesWidget is used
        for identifier, code_model in self.get_models():
            # Skip if the code is not displayed, convenient for the plugin developer
            if not code_model.is_ready:
                continue
            if not issubclass(
                code_model.code_widget_class, QEAppComputationalResourcesWidget
            ):
                yield (
                    f"Error: hi, plugin developer, please use the QEAppComputationalResourcesWidget from aiidalab_qe.common.widgets for code {identifier}."
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
            Reference time of a single scf_cycle. Default is 129.6.
        `tmax` : `float`, optional
            Maximum time limit. Default is 12 hours.
        `scf_cycles` : `int`, optional
            Reference number of SCF cycles in a relaxation. Default is 5.

        The default values n0, v0, num_cpus0, t0 and scf_cycles are taken from the simulation of SiO2 bulk
        example structure present in the app, following a moderate protocol. We then used the formula

        num_cpus = num_cpus0 * (n/n0)^3 * (v/v0)^(3/2) * (scf_cycles * t0)/tmax

        assuming that the number of CPUs scales with the number of atoms as power of 3, the volume of the system as power of 3/2,
        the number of SCF cycles and the time of a single SCF cycle. The power dependence was then adjusted to match the
        other reference calculations done on bulk SiO2, Silicon and Gold, using different number of cpus.
        NOTE: this is a very rough estimate and should be used as a guideline only.

        Returns
        -------
        `int`
            The estimated minimum number of CPUs required.
        """
        import numpy as np

        return int(
            np.ceil(
                num_cpus0
                * (n / n0) ** 3
                * (v / v0) ** (3 / 2)
                * (scf_cycles * t0)
                / tmax
            )
        )
