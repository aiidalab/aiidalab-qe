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
        return "relax" in self.properties

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
                <strong>PK:</strong> {structure.pk}<br>
                <strong>Label:</strong> {structure.label}<br>
                <strong>Description:</strong> {structure.description}<br>
                <strong>Number of atoms:</strong> {len(structure.sites)}<br>
                <strong>Creation time:</strong> {formatted} ({relative})<br>
            </div>
        """

    def _get_atom_table_data(self):
        """Build table data; if 'fixed_atoms' is present in structure attributes,
        add a 'Free x,y,z' column showing '✓' for free (1) and 'x' for fixed (0)."""
        # Try to get fixed_atoms from AiiDA StructureData attributes
        fixed_atoms = None
        try:
            fixed_atoms = self.structure.base.attributes.all['fixed_atoms']
        except KeyError:
            fixed_atoms = None

        ase_atoms = self.structure.get_ase()

        # Header
        data = [["Atom index", "Chemical symbol", "Tag", "x (Å)", "y (Å)", "z (Å)"]]
        if fixed_atoms is not None:
            data[0].append("Free x,y,z")

        positions = ase_atoms.positions
        chemical_symbols = ase_atoms.get_chemical_symbols()
        tags = ase_atoms.get_tags()

        def fmt_free(mask):
            """mask: tuple/list of three 0/1; 0=free -> 'X', 1=fixed -> ' '."""
            try:
                x, y, z = mask
            except Exception:
                x = y = z = 0
            # If your UI collapses spaces, replace ' ' with '·' or '\u00A0' (NBSP).
            return f"({'x' if x == 0 else '✓'} {'x' if y == 0 else '✓'} {'x' if z == 0 else '✓'})"

        for idx, (symbol, tag, pos) in enumerate(zip(chemical_symbols, tags, positions), start=1):
            formatted_position = [f"{coord:.2f}" for coord in pos]
            row = [idx, symbol, tag, *formatted_position]

            if fixed_atoms is not None:
                mask = fixed_atoms[idx - 1] if idx - 1 < len(fixed_atoms) else (0, 0, 0)
                row.append(fmt_free(mask))

            data.append(row)

        return data


    def get_model_state(self):
        return {
            "selected_view": self.selected_view,
        }

    def set_model_state(self, parameters):
        self.selected_view = parameters.get("selected_view", "initial")
