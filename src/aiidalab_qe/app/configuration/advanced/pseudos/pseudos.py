from __future__ import annotations

from copy import deepcopy

import ipywidgets as ipw

from aiida import orm
from aiidalab_qe.common.panel import ConfigurationSettingsPanel
from aiidalab_qe.common.widgets import HBoxWithUnits
from aiidalab_qe.utils import generate_alert
from aiidalab_widgets_base.utils import StatusHTML

from .model import PseudosConfigurationSettingsModel
from .uploader import PseudoPotentialUploader, PseudoPotentialUploaderModel


class PseudosConfigurationSettingsPanel(
    ConfigurationSettingsPanel[PseudosConfigurationSettingsModel],
):
    def __init__(self, model: PseudosConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_protocol_change,
            "protocol",
        )
        self._model.observe(
            self._on_spin_orbit_change,
            "spin_orbit",
        )
        self._model.observe(
            self._on_family_parameters_change,
            ["library", "functional"],
        )
        self._model.observe(
            self._on_functionals_change,
            "functionals",
        )
        self._model.observe(
            self._on_family_change,
            "family",
        )
        self._model.observe(
            self._on_dictionary_change,
            "dictionary",
        )
        self._model.observe(
            self._on_cutoffs_change,
            "cutoffs",
        )
        self._model.observe(
            self._on_show_upload_warning_change,
            "show_upload_warning",
        )

    def render(self):
        if self.rendered:
            return

        self.family_header = ipw.HTML()
        ipw.dlink(
            (self._model, "family_header"),
            (self.family_header, "value"),
        )

        self.family_help = ipw.HTML()
        ipw.dlink(
            (self._model, "family_help_message"),
            (self.family_help, "value"),
        )

        self.functional = ipw.ToggleButtons()
        ipw.dlink(
            (self._model, "functional_options"),
            (self.functional, "options"),
        )
        ipw.link(
            (self._model, "functional"),
            (self.functional, "value"),
        )

        self.library = ipw.ToggleButtons(style={"button_width": "fit-content"})
        ipw.dlink(
            (self._model, "library_options"),
            (self.library, "options"),
        )
        ipw.link(
            (self._model, "library"),
            (self.library, "value"),
        )

        self.setter_widget = ipw.VBox()

        self.ecutwfc = ipw.FloatText(
            description="Wavefunction",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "ecutwfc"),
            (self.ecutwfc, "value"),
        )
        self.ecutrho = ipw.FloatText(
            description="Charge density",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "ecutrho"),
            (self.ecutrho, "value"),
        )

        self._status_message = StatusHTML(clear_after=20)
        ipw.link(
            (self._model, "status_message"),
            (self._status_message, "message"),
        )

        self._warning_message = ipw.HTML(
            value=generate_alert(
                alert_type="warning",
                message="""
                    <h4><b>⚠️ Pseudopotential upload detected ⚠️</b></h4>
                    Adjust plane-wave cutoff energies below to ensure convergence
                    <br>
                    To reset to protocol defaults, select both an exchange-correlation
                    functional and a pseudopotential family above.
                    <br>
                    Switching the protocol in <b>Basic settings</b> will also reset the
                    pseudopotentials to the defaults of the selected protocol.
                """,
                class_="text-center",
                style_="line-height: 2;",
            ),
            layout=ipw.Layout(
                display="block" if self._model.show_upload_warning else "none",
            ),
        )

        self.children = [
            ipw.HTML("""
                <div class="pseudo-text">
                    The exchange-correlation functional and pseudopotential library are
                    set by the <b>protocol</b> configured in the <b>Basic settings</b>
                    tab.
                    <br>
                    Here you can override the defaults if desired.
                </div>
            """),
            ipw.HTML("<h4>Exchange-correlation functional</h4>"),
            ipw.HTML("""
                <div class="pseudo-text">
                    The exchange-correlation energy is described using the selected
                    functional.
                    <br>
                    We currently provide support for two well-established generalized
                    gradient approximation (GGA) functionals: PBE and PBEsol.
                </div>
            """),
            self.functional,
            self.family_header,
            self.family_help,
            self.library,
            ipw.HTML("<h4>Pseudopotentials</h4>"),
            self._warning_message,
            ipw.HTML("""
                <div class="pseudo-text">
                    The pseudopotential for each kind of atom in the structure can be
                    custom set.
                    <br>
                    The default pseudopotential and cutoffs are taken from the
                    pseudopotential family.
                    <br>
                    Pseudopotential information is shown (if provided) in the following order:
                    <ul>
                        <li>
                            Recommended wavefunction cutoff energy (Ry)
                        </li>
                        <li>
                            Recommended charge density cutoff energy (Ry)
                        </li>
                        <li>
                            Exchange-correlation functional
                        </li>
                        <li>
                            Relativistic treatment
                        </li>
                    </ul>
                </div>
            """),
            self.setter_widget,
            ipw.HTML("<h4>Cutoffs</h4>"),
            ipw.HTML("""
                <div style="line-height: 1.4;">
                    The
                    <a
                        href="https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm312"
                        target="_blank"
                    >
                        default cutoffs
                    </a> used for the calculation are the maximum of the default cutoffs
                    from all pseudopotentials.
                    <br>
                    You can override them here.
                </div>
            """),
            HBoxWithUnits(self.ecutwfc, "Ry"),
            HBoxWithUnits(self.ecutrho, "Ry"),
            self._status_message,
        ]

        self.rendered = True

        self.refresh(specific="widgets")

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")

    def _on_spin_orbit_change(self, _):
        self._model.update_library_options()

    def _on_family_parameters_change(self, _):
        self._model.update_family()
        self._model.update_family_header()
        self._model.update_functionals()
        self._model.update_blockers()

    def _on_functionals_change(self, _):
        self._model.functional = (
            self._model.functionals[0]  # type: ignore
            if len(set(self._model.functionals)) == 1
            else None
        )
        self._model.update_blockers()

    def _on_family_change(self, _):
        self._model.show_upload_warning = not self._model.family
        self._model.update_dictionary()

    def _on_dictionary_change(self, _):
        self._model.update_cutoffs()

    def _on_cutoffs_change(self, change):
        cutoffs = change["new"]  # [[ecutwfc...], [ecutrho...]]
        self._model.ecutwfc = max(cutoffs[0])
        self._model.ecutrho = max(cutoffs[1])

    def _on_show_upload_warning_change(self, change):
        if not self.rendered:
            return
        self._warning_message.layout.display = "block" if change["new"] else "none"

    def update(self, specific=""):
        if self._model.updated:
            return
        self._show_loading()
        if not self._model.locked or (specific and specific != "widgets"):
            self._model.update(specific)
        self._build_setter_widgets()
        self._model.update_library_options()
        self._model.update_family_header()
        self._model.updated = True

    def _show_loading(self):
        if self.rendered:
            self.setter_widget.children = [self.loading_message]

    def _build_setter_widgets(self):
        if not self.rendered:
            return

        children = []

        kinds = self._model.input_structure.kinds if self._model.has_structure else []

        for index, kind in enumerate(kinds):
            uploader_model = PseudoPotentialUploaderModel(kind.name, kind.symbol)
            uploader = PseudoPotentialUploader(model=uploader_model)

            def on_default_pseudo(
                _=None,
                model: PseudoPotentialUploaderModel = uploader_model,
                kind_name=kind.name,
                index=index,
            ):
                if not self._model.dictionary:
                    return
                try:
                    pseudo = orm.load_node(self._model.dictionary.get(kind_name))
                except Exception:
                    pseudo = None
                model.pseudo = pseudo
                model.cutoffs = (
                    [
                        self._model.cutoffs[0][index],  # type: ignore
                        self._model.cutoffs[1][index],  # type: ignore
                    ]
                    if len(self._model.cutoffs[0]) > index  # type: ignore
                    else [0.0, 0.0]
                )
                model.update_pseudo_info()
                model.uploaded = False

            self._model.observe(
                on_default_pseudo,
                "dictionary",
            )

            def on_pseudo_upload(
                change,
                model: PseudoPotentialUploaderModel = uploader_model,
                kind_name=kind.name,
                index=index,
            ):
                if not (change["new"] and model.pseudo):
                    return

                self._model.family = None
                self._model.library = None
                self._model.functional = None

                self._model.dictionary = {
                    **self._model.dictionary,
                    kind_name: model.pseudo.uuid,
                }

                functional = model.pseudo.base.extras.get("functional", None)
                functionals = [*self._model.functionals]
                functionals[index] = functional
                # The following double-setting is done to force the blockers check,
                # which now also bail early if the functionals are empty.
                self._model.functionals = []
                self._model.functionals = functionals

                cutoffs: list = deepcopy(self._model.cutoffs)  # type: ignore
                cutoffs[0][index] = model.cutoffs[0]  # type: ignore
                cutoffs[1][index] = model.cutoffs[1]  # type: ignore
                self._model.cutoffs = cutoffs

            uploader_model.observe(
                on_pseudo_upload,
                "uploaded",
            )

            uploader.render()

            on_default_pseudo()

            self._links.extend(uploader.links)

            children.append(uploader)

        self.setter_widget.children = children
