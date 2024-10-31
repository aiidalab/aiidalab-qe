from __future__ import annotations

import os
import typing as t
from copy import deepcopy

import traitlets as tl

from aiida import orm
from aiida.engine import ProcessBuilderNamespace
from aiida.engine import submit as aiida_submit
from aiida.orm.utils.serialize import serialize
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget
from aiidalab_qe.workflows import QeAppWorkChain

from .code import CodeModel, CodesDict

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmissionModel(tl.HasTraits):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    input_parameters = tl.Dict()

    process = tl.Instance(orm.WorkChainNode, allow_none=True)
    process_label = tl.Unicode("")
    process_description = tl.Unicode("")

    submission_blocker_messages = tl.Unicode("")
    submission_warning_messages = tl.Unicode("")

    installing_qe = tl.Bool(False)
    installing_sssp = tl.Bool(False)
    qe_installed = tl.Bool(allow_none=True)
    sssp_installed = tl.Bool(allow_none=True)

    codes = tl.Dict(
        key_trait=tl.Unicode(),  # plugin identifier
        value_trait=tl.Dict(  # plugin codes
            key_trait=tl.Unicode(),  # code name
            value_trait=tl.Instance(CodeModel),  # code metadata
        ),
        default_value={},
    )

    internal_submission_blockers = tl.List(tl.Unicode())
    external_submission_blockers = tl.List(tl.Unicode())

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

        # Used by the code-setup thread to fetch code options
        # This is necessary to avoid passing the User object
        # between session in separate threads.
        self._default_user_email = orm.User.collection.get_default().email

    @property
    def is_blocked(self):
        return any(
            [
                *self.internal_submission_blockers,
                *self.external_submission_blockers,
            ]
        )

    def submit(self):
        parameters = self.get_model_state()
        builder = self._create_builder(parameters)

        with self.hold_trait_notifications():
            process = aiida_submit(builder)

            process.label = self.process_label
            process.description = self.process_description
            # since AiiDA data node may exist in the ui_parameters,
            # we serialize it to yaml
            process.base.extras.set("ui_parameters", serialize(parameters))
            # store the workchain name in extras, this will help to filter the workchain in the future
            process.base.extras.set("workchain", parameters["workchain"])  # type: ignore
            process.base.extras.set(
                "structure",
                self.input_structure.get_formula(),
            )
            self.process = process

    def check_resources(self):
        pw_code = self.get_code("dft", "pw")

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

        self.submission_warning_messages = (
            ""
            if (on_localhost and num_cpus / localhost_cpus) <= 0.8
            and (not large_system or estimated_CPUs <= num_cpus)
            else self._ALERT_MESSAGE.format(
                alert_class="warning",
                message=alert_message,
            )
        )

    def refresh_codes(self):
        for _, code_model in self.get_code_models(flat=True):
            code_model.update(self._default_user_email)  # type: ignore

    def update_active_codes(self):
        for name, code_model in self.get_code_models(flat=True):
            if name != "pw":
                code_model.deactivate()
        properties = self._get_properties()
        for identifier, code_models in self.get_code_models():
            if identifier in properties:
                for code_model in code_models.values():
                    code_model.activate()

    def update_process_label(self):
        if not self.input_structure:
            self.process_label = ""
            return
        structure_label = (
            self.input_structure.label
            if len(self.input_structure.label) > 0
            else self.input_structure.get_formula()
        )
        workchain_data = self.input_parameters.get(
            "workchain",
            {"properties": []},
        )
        properties = [p for p in workchain_data["properties"] if p != "relax"]
        relax_type = workchain_data.get("relax_type", "none")
        relax_info = "unrelaxed"
        if relax_type != "none":
            relax_info = (
                "relax: atoms+cell" if "cell" in relax_type else "relax: atoms only"
            )
        protocol_and_magnetic_info = f"{workchain_data['protocol']} protocol"
        if workchain_data["spin_type"] != "none":
            protocol_and_magnetic_info += ", magnetic"
        properties_info = f"→ {', '.join(properties)}" if properties else ""
        label = f"{structure_label} [{relax_info}, {protocol_and_magnetic_info}] {properties_info}".strip()
        self.process_label = label

    def update_submission_blockers(self):
        self.internal_submission_blockers = list(self._check_submission_blockers())

    def update_submission_blocker_message(self):
        blockers = self.internal_submission_blockers + self.external_submission_blockers
        if any(blockers):
            fmt_list = "\n".join(f"<li>{item}</li>" for item in sorted(blockers))
            self.submission_blocker_messages = f"""
                <div class="alert alert-info">
                    <b>The submission is blocked due to the following reason(s):</b>
                    <ul>
                        {fmt_list}
                    </ul>
                </div>
            """
        else:
            self.submission_blocker_messages = ""

    def get_model_state(self) -> dict[str, dict[str, dict]]:
        parameters: dict = deepcopy(self.input_parameters)  # type: ignore
        parameters["codes"] = self.get_selected_codes()
        return parameters

    def set_model_state(self, parameters):
        if "resources" in parameters:
            parameters["codes"] = {
                key: {"code": value} for key, value in parameters["codes"].items()
            }
            parameters["codes"]["pw"]["nodes"] = parameters["resources"]["num_machines"]
            parameters["codes"]["pw"]["cpus"] = parameters["resources"][
                "num_mpiprocs_per_machine"
            ]
            parameters["codes"]["pw"]["parallelization"] = {
                "npool": parameters["resources"]["npools"]
            }
        self.set_selected_codes(parameters["codes"])
        if self.process:
            self.process_label = self.process.label
            self.process_description = self.process.description

    def add_code(self, identifier, name, code):
        code.name = name
        if identifier not in self.codes:
            self.codes[identifier] = {}  # type: ignore
        self.codes[identifier][name] = code  # type: ignore

    def get_code(self, identifier, name) -> CodeModel | None:
        if identifier in self.codes and name in self.codes[identifier]:  # type: ignore
            return self.codes[identifier][name]  # type: ignore

    def get_code_models(
        self,
        flat=False,
    ) -> t.Iterator[tuple[str, CodesDict | CodeModel]]:
        if flat:
            for codes in self.codes.values():
                yield from codes.items()
        else:
            yield from self.codes.items()

    def get_selected_codes(self) -> dict[str, dict]:
        return {
            name: code_model.get_model_state()
            for name, code_model in self.get_code_models(flat=True)
            if code_model.is_ready
        }

    def set_selected_codes(self, code_data=DEFAULT["codes"]):
        with self.hold_trait_notifications():
            for name, code_model in self.get_code_models(flat=True):
                if name in code_data:
                    code_model.set_model_state(code_data[name])

    def reset(self):
        with self.hold_trait_notifications():
            self.input_structure = None
            self.input_parameters = {}
            self.process = None

    def _get_properties(self) -> list[str]:
        return self.input_parameters.get("workchain", {}).get("properties", [])

    def _create_builder(self, parameters) -> ProcessBuilderNamespace:
        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            parameters=deepcopy(parameters),  # TODO why deepcopy again?
        )

        codes = parameters["codes"]

        builder.relax.base.pw.metadata.options.resources = {
            "num_machines": codes.get("pw")["nodes"],
            "num_mpiprocs_per_machine": codes.get("pw")["ntasks_per_node"],
            "num_cores_per_mpiproc": codes.get("pw")["cpus_per_task"],
        }
        mws = codes.get("pw")["max_wallclock_seconds"]
        builder.relax.base.pw.metadata.options["max_wallclock_seconds"] = mws
        parallelization = codes["pw"]["parallelization"]
        builder.relax.base.pw.parallelization = orm.Dict(dict=parallelization)

        return builder

    def _check_submission_blockers(self):
        # Do not submit while any of the background setup processes are running.
        if self.installing_qe or self.installing_sssp:
            yield "Background setup processes must finish."

        # SSSP library not installed
        if not self.sssp_installed:
            yield "The SSSP library is not installed."

        # No pw code selected (this is ignored while the setup process is running).
        pw_code = self.get_code(identifier="dft", name="pw")
        if pw_code and not pw_code.selected and not self.installing_qe:
            yield ("No pw code selected")

        # code related to the selected property is not installed
        properties = self._get_properties()
        message = "Calculating the {property} property requires code {code} to be set."
        for identifier, codes in self.get_code_models():
            if identifier in properties:
                for code in codes.values():
                    if not code.is_ready:
                        yield message.format(property=identifier, code=code.description)

        # check if the QEAppComputationalResourcesWidget is used
        for name, code in self.get_code_models(flat=True):
            # skip if the code is not displayed, convenient for the plugin developer
            if not code.is_ready:
                continue
            if not issubclass(
                code.code_widget_class, QEAppComputationalResourcesWidget
            ):
                yield (
                    f"Error: hi, plugin developer, please use the QEAppComputationalResourcesWidget from aiidalab_qe.common.widgets for code {name}."
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
