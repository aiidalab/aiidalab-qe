from __future__ import annotations

import numpy as np
import traitlets as tl

from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.common.mixins import HasInputStructure

from ..subsettings import AdvancedCalculationSubSettingsModel


class ConvergenceConfigurationSettingsModel(
    AdvancedCalculationSubSettingsModel,
    HasInputStructure,
):
    title = "Convergence"
    identifier = "convergence"

    dependencies = [
        "input_structure",
        "protocol",
    ]

    protocol = tl.Unicode()

    help_message = tl.Unicode()
    scf_conv_thr = tl.Float(0.0)
    scf_conv_thr_step = tl.Float(1e-10)
    etot_conv_thr = tl.Float(0.0)
    etot_conv_thr_step = tl.Float(1e-5)
    forc_conv_thr = tl.Float(0.0)
    forc_conv_thr_step = tl.Float(1e-4)
    electron_maxstep = tl.Int(80)
    optimization_maxsteps = tl.Int(50)
    mixing_mode_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "plain",
            "TF",
            "local-TF",
        ],
    )
    mixing_mode = tl.Unicode("plain")
    mixing_beta = tl.Float(0.4)

    def update(self, specific=""):
        if specific == "structure":
            self._update_help_message()
        elif specific == "protocol":
            self._update_thresholds()

    def reset(self):
        with self.hold_trait_notifications():
            self.scf_conv_thr = self._get_default("scf_conv_thr")
            self.scf_conv_thr_step = self._get_default("scf_conv_thr_step")
            self.etot_conv_thr = self._get_default("etot_conv_thr")
            self.etot_conv_thr_step = self._get_default("etot_conv_thr_step")
            self.forc_conv_thr = self._get_default("forc_conv_thr")
            self.forc_conv_thr_step = self._get_default("forc_conv_thr_step")
            self.electron_maxstep = self._get_default("electron_maxstep")
            self.optimization_maxsteps = self._get_default("optimization_maxsteps")
            self.mixing_mode = self._get_default("mixing_mode")
            self.mixing_beta = self._get_default("mixing_beta")

    def _update_help_message(self):
        if not self.has_structure:
            self.help_message = "No structure available."
            return
        num_atoms = len(self.input_structure.sites)
        self.help_message = f"""
            <div style="line-height: 1.4; margin-bottom: 5px;">
                Setting the energy threshold for the self-consistent field (SCF)
                and energy and force thresholds for ionic convergence ensures
                calculation accuracy and stability. Lower values increase the
                accuracy but also the computational cost. The default values set by
                the <b>protocol</b> are usually a good starting point. For energy
                thresholds, the actual value used in the calculation (shown below
                widget) is given as:
                <code>threshold x num_atoms</code>
                (<code>num_atoms = {num_atoms}</code>)
            </div>
        """

    def _update_thresholds(self):
        parameters = PwBaseWorkChain.get_protocol_inputs(self.protocol)

        etot_value = parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        self._set_value_and_step("etot_conv_thr", etot_value)
        self.etot_conv_thr = self._defaults["etot_conv_thr"]
        self.etot_conv_thr_step = self._defaults["etot_conv_thr_step"]

        scf_value = parameters["meta_parameters"]["conv_thr_per_atom"]
        self._set_value_and_step("scf_conv_thr", scf_value)
        self.scf_conv_thr = self._defaults["scf_conv_thr"]
        self.scf_conv_thr_step = self._defaults["scf_conv_thr_step"]

        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        self._set_value_and_step("forc_conv_thr", forc_value)
        self.forc_conv_thr = self._defaults["forc_conv_thr"]
        self.forc_conv_thr_step = self._defaults["forc_conv_thr_step"]

    def _set_value_and_step(self, attribute, value):
        self._defaults[attribute] = value
        if value != 0:
            order_of_magnitude = np.floor(np.log10(abs(value)))
            step = 10 ** (order_of_magnitude - 1)
        else:
            step = 0.1
        self._defaults[f"{attribute}_step"] = step

    def _get_default(self, trait):
        return self._defaults.get(trait, self.traits()[trait].default_value)
