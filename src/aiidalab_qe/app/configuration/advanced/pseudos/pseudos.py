from __future__ import annotations

import typing as t

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
            "structure_uuid",
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
            self._on_functional_change,
            "functional",
        )
        self._model.observe(
            self._on_library_change,
            "library",
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
            self._on_ecut_change,
            ["ecutwfc", "ecutrho"],
        )

        self._uploader_dictionary_observers: list[t.Callable] = []

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

        self.pseudos_list = ipw.VBox()

        self.ecutwfc = ipw.BoundedFloatText(
            description="Wavefunction",
            min=0.0,
            max=1e4,
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "ecutwfc"),
            (self.ecutwfc, "value"),
        )
        self.ecutrho = ipw.BoundedFloatText(
            description="Charge density",
            min=0.0,
            max=1e4,
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
        ipw.dlink(
            (self._model, "show_upload_warning"),
            (self._warning_message.layout, "display"),
            lambda show: "block" if show else "none",
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
            self.pseudos_list,
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
        self._unsubscribe_uploader_observations()
        self.refresh(specific="structure")

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")

    def _on_spin_orbit_change(self, _):
        self._model.update_library_options()

    def _on_functional_change(self, _):
        self._model.update_family()

    def _on_library_change(self, _):
        self._model.update_family_header()
        self._model.update_family()

    def _on_family_change(self, _):
        self._toggle_upload_warning()
        self._model.update_dictionary()

    def _on_dictionary_change(self, _):
        self._model.update_functional()
        self._model.update_cutoffs()
        self._model.update_blockers()

    def _on_cutoffs_change(self, change):
        cutoffs = change["new"]  # [[ecutwfc...], [ecutrho...]]
        self._model.ecutwfc = max(cutoffs[0])
        self._model.ecutrho = max(cutoffs[1])

    def _on_ecut_change(self, _):
        self._model.update_blockers()

    def _update_ui(self):
        self._build_pseudos_list()
        self._model.update_family_header()

    def _unsubscribe_uploader_observations(self):
        for observer in self._uploader_dictionary_observers:
            self._model.unobserve(observer, "dictionary")
        self._uploader_dictionary_observers.clear()

    def _toggle_upload_warning(self):
        self._model.show_upload_warning = not self._model.family

    def _build_pseudos_list(self):
        if not self.rendered:
            return

        self.pseudos_list.children = [self.loading_message]

        children = []

        kinds = self._model.input_structure.kinds if self._model.has_structure else []

        for index, kind in enumerate(kinds):
            uploader_model = PseudoPotentialUploaderModel(kind.name, kind.symbol)

            def on_pseudo_upload(
                _=None,
                model=uploader_model,
                kind_name=kind.name,
            ):
                if not (model.uploaded and model.pseudo):
                    return
                with self._model.hold_trait_notifications():
                    self._model.dictionary = {
                        **self._model.dictionary,
                        kind_name: model.pseudo.uuid,
                    }
                    self._model.library = None
                    self._model.family = None

            uploader_model.observe(
                on_pseudo_upload,
                "uploaded",
            )

            def synchronize_pseudo(
                _=None,
                model=uploader_model,
                kind_name=kind.name,
                index=index,
            ):
                if not self._model.dictionary:
                    return
                try:
                    pseudo = orm.load_node(self._model.dictionary.get(kind_name))
                except Exception:
                    pseudo = None
                model.uploaded = model.pseudo is pseudo
                model.pseudo = pseudo
                model.cutoffs = self._model.get_cutoffs_by_index(index)
                model.update_pseudo_info()

            self._model.observe(
                synchronize_pseudo,
                "dictionary",
            )
            self._uploader_dictionary_observers.append(synchronize_pseudo)

            uploader = PseudoPotentialUploader(model=uploader_model)
            uploader.render()

            self._links.extend(uploader.links)

            children.append(uploader)

            synchronize_pseudo()

        self.pseudos_list.children = children
