from __future__ import annotations
from pydoc import html

from aiida_quantumespresso.tools import data
from matplotlib import legend
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

    # def _get_atom_table_data(self):
    #     """Build table data; if 'fixed_atoms' is present in structure attributes,
    #     add a 'Free x,y,z' column showing '✓' for free (1) and 'x' for fixed (0)."""
    #     # Try to get fixed_atoms from AiiDA StructureData attributes
    #     fixed_atoms = self.structure.base.attributes.all.get('fixed_atoms', None) 

    #     ase_atoms = self.structure.get_ase()

    #     # Header
    #     data = [["Atom index", "Chemical symbol", "Tag", "x (Å)", "y (Å)", "z (Å)"]]
    #     if fixed_atoms is not None:
    #         data[0].append("Free x,y,z")

    #     positions = ase_atoms.positions
    #     chemical_symbols = ase_atoms.get_chemical_symbols()
    #     tags = ase_atoms.get_tags()

    #     def fmt_free(mask):
    #         """mask: tuple/list of three 0/1; 0=free -> 'X', 1=fixed -> ' '."""
    #         try:
    #             x, y, z = mask
    #         except Exception:
    #             x = y = z = 0
    #         # If your UI collapses spaces, replace ' ' with '·' or '\u00A0' (NBSP).
    #         return f"({'x' if x == 0 else '✓'} {'x' if y == 0 else '✓'} {'x' if z == 0 else '✓'})"

    #     for idx, (symbol, tag, pos) in enumerate(zip(chemical_symbols, tags, positions), start=1):
    #         formatted_position = [f"{coord:.2f}" for coord in pos]
    #         row = [idx, symbol, tag, *formatted_position]

    #         if fixed_atoms is not None:
    #             mask = fixed_atoms[idx - 1] if idx - 1 < len(fixed_atoms) else (0, 0, 0)
    #             row.append(fmt_free(mask))

    #         data.append(row)

    #     return data
    def _get_atom_table_data(self):
        """Build table data.

        - If 'fixed_atoms' is present in AiiDA StructureData attributes (or ASE arrays),
        add a 'Free x,y,z' column showing '✓' for free (1) and 'x' for fixed (0).
        - If constraints exist (CONSTRAINTS block), add a 'Constr.' column with
        star labels (*1, *2, …) for atoms involved in each constraint.
        """
        # --- source data ---
        attrs_all = getattr(getattr(getattr(self.structure, "base", None), "attributes", None), "all", {}) or {}
        fixed_atoms = attrs_all.get("fixed_atoms", None)

        ase_atoms = self.structure.get_ase()

        # Try also to read fixed mask from ASE arrays if attributes missing
        if fixed_atoms is None and "fixed_atoms" in getattr(ase_atoms, "arrays", {}):
            fixed_atoms = ase_atoms.arrays["fixed_atoms"].tolist()

        # constraints: prefer attributes, else ASE info
        constraints_info = attrs_all.get("CONSTRAINTS")
        if constraints_info is None:
            constraints_info = ase_atoms.info.get("CONSTRAINTS", None) if hasattr(ase_atoms, "info") else None
        # normalize constraints_info
        if not isinstance(constraints_info, dict):
            constraints_info = {"number": 0, "tolerance": "1e-6", "list": []}
        constraints_list = constraints_info.get("list", []) or []

        # Build star marks per atom and a legend (legend not returned here)
        def build_marks(n_atoms, entries):
            marks = {i: [] for i in range(n_atoms)}
            legend = []
            for k, entry in enumerate(entries, start=1):
                legend.append(f"*{k} {entry}")
                parts = str(entry).strip().split()
                if len(parts) >= 4 and parts[0] == "distance":
                    try:
                        i1 = int(parts[1]) - 1
                        i2 = int(parts[2]) - 1
                    except Exception:
                        continue
                    label = f"*{k}"
                    if 0 <= i1 < n_atoms:
                        marks[i1].append(label)
                    if 0 <= i2 < n_atoms:
                        marks[i2].append(label)
            return marks, legend

        marks_map, legend = build_marks(len(ase_atoms), constraints_list)

        # --- header ---
        data = [["Atom index", "Chemical symbol", "Tag", "x (Å)", "y (Å)", "z (Å)"]]
        if fixed_atoms is not None:
            data[0].append("Free x,y,z")
        # add constraints column if any constraints exist
        if constraints_list:
            data[0].append("Constr.")

        positions = ase_atoms.positions
        chemical_symbols = ase_atoms.get_chemical_symbols()
        tags = ase_atoms.get_tags()

        def fmt_free(mask):
            """mask: 3-tuple/list of 0/1; 1=free → '✓', 0=fixed → 'x'."""
            try:
                x, y, z = mask
            except Exception:
                x = y = z = 1  # default free if malformed
            return f"({'✓' if int(x)==1 else 'x'} {'✓' if int(y)==1 else 'x'} {'✓' if int(z)==1 else 'x'})"

        for idx, (symbol, tag, pos) in enumerate(zip(chemical_symbols, tags, positions), start=1):
            formatted_position = [f"{coord:.2f}" for coord in pos]
            row = [idx, symbol, tag, *formatted_position]

            if fixed_atoms is not None:
                # robustly fetch mask; default to all free if missing
                mask = (1, 1, 1)
                try:
                    if 0 <= (idx - 1) < len(fixed_atoms):
                        mask = fixed_atoms[idx - 1]
                except Exception:
                    pass
                row.append(fmt_free(mask))

            if constraints_list:
                stars = " ".join(marks_map.get(idx - 1, []))
                row.append(stars)

            data.append(row)
        # temporary quick fix to show teh legend of constraints
        if legend:
            ncols = len(data[0])
            legend_text = "; ".join(legend)

            injected = (
                # close the last normal row
                f'</td></tr>'
                # add a non-interactive full-width legend row
                f'<tr class="constraints-legend-row" '
                f'style="pointer-events:none; user-select:none;">'
                f'<td colspan="{ncols}" style="text-align:right;">{html.escape(legend_text)}</td>'
                f'</tr>'
                # open a dummy row/cell so the renderer’s trailing </td></tr> stays balanced
                f'<tr><td>'
            )
            data.append([injected])



        return data


    def get_model_state(self):
        return {
            "selected_view": self.selected_view,
        }

    def set_model_state(self, parameters):
        self.selected_view = parameters.get("selected_view", "initial")
