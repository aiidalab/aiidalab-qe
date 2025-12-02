from __future__ import annotations

import traitlets as tl
from ase.formula import Formula

from aiida import orm
from aiidalab_qe.common.panel import ResultsModel
from aiidalab_qe.common.time import format_time, relative_time


class StructureResultsModel(ResultsModel):
    title = "Structure"
    identifier = "structure"

    structure = tl.Instance(orm.StructureData, allow_none=True)
    selected_view = tl.Unicode("initial")
    header = tl.Unicode()
    sub_header = tl.Unicode()
    source = tl.Instance(orm.utils.managers.NodeLinksManager, allow_none=True)
    info = tl.Unicode()
    table_data = tl.List(tl.List())

    _this_process_label = "PwRelaxWorkChain"

    HEADER_TEMPLATE = "<h2 style='margin: 0;'>{content}</h2>"
    _SUB_HEADER_TEMPLATE = "<h4 style='margin: 0; font-size: 18px'>{content}</h4>"

    @property
    def include(self):
        return True

    @property
    def is_relaxed(self):
        if "relax" not in self.properties:
            return False
        parameters = self.inputs.relax.base.pw.parameters.get_dict()
        return "relax" in parameters["CONTROL"]["calculation"]

    def update(self):
        super().update()
        with self.hold_trait_notifications():
            if not self.is_relaxed or self.selected_view == "initial":
                self.sub_header = self._SUB_HEADER_TEMPLATE.format(content="Initial")
                self.source = self.inputs
            else:
                self.sub_header = self._SUB_HEADER_TEMPLATE.format(content="Relaxed")
                self.source = self.outputs
            self.structure = self._get_structure()
            if self.structure:
                self.info = self._get_structure_info()
                self.table_data = self._get_atom_table_data()
                if not self.header:
                    formatted = Formula(self.structure.get_formula()).format("html")
                    self.header = self.HEADER_TEMPLATE.format(content=formatted)

    def toggle_selected_view(self):
        self.selected_view = "relaxed" if self.selected_view == "initial" else "initial"

    def _get_structure(self) -> orm.StructureData | None:
        try:
            return self.source.structure if self.source else None
        except AttributeError:
            # If source is outputs but job failed, there may not be a structure
            return None

    def _get_structure_info(self):
        structure = self.structure
        formatted = format_time(structure.ctime)
        relative = relative_time(structure.ctime)
        return f"""
            <div style='line-height: 1.4;'>
                <b>PK:</b> {structure.pk}<br>
                <b>Label:</b> {structure.label}<br>
                <b>Description:</b> {structure.description}<br>
                <b>Number of atoms:</b> {len(structure.sites)}<br>
                <b>Creation time:</b> {formatted} ({relative})<br>
            </div>
        """

    def _get_atom_table_data(self):
        structure = self.structure.get_ase()
        data = [
            [
                "Atom index",
                "Chemical symbol",
                "Tag",
                "x (Å)",
                "y (Å)",
                "z (Å)",
            ]
        ]
        positions = structure.positions
        chemical_symbols = structure.get_chemical_symbols()
        tags = structure.get_tags()

        for index, (symbol, tag, position) in enumerate(
            zip(chemical_symbols, tags, positions), start=1
        ):
            formatted_position = [f"{coord:.2f}" for coord in position]
            data.append([index, symbol, tag, *formatted_position])

        return data

    def get_model_state(self):
        return {
            "selected_view": self.selected_view,
        }

    def set_model_state(self, parameters):
        self.selected_view = parameters.get("selected_view", "initial")
