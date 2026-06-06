from __future__ import annotations

import typing as t

import traitlets as tl

from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

if t.TYPE_CHECKING:
    from aiidalab_qe.common.mvc import Model

StructureType = t.Union[orm.StructureData, HubbardStructureData]


class HasInputStructure(tl.HasTraits):
    structure_uuid = tl.Unicode(None, allow_none=True)

    @property
    def input_structure(self) -> StructureType | None:
        if not self.structure_uuid:
            return None
        try:
            return t.cast(StructureType, orm.load_node(self.structure_uuid))
        except NotExistent:
            return None

    @property
    def has_structure(self):
        return self.input_structure is not None

    @property
    def has_pbc(self):
        return not self.has_structure or any(self.input_structure.pbc)

    @property
    def has_tags(self):
        return self.has_structure and any(
            not kind_name.isalpha()
            for kind_name in self.input_structure.get_kind_names()
        )

    def get_num_atoms(self) -> int:
        return len(self.input_structure.sites) if self.has_structure else 0


M = t.TypeVar("M", bound="Model")


class HasModels(t.Generic[M]):
    def __init__(self):
        self._models: dict[str, M] = {}

    def has_model(self, identifier: str):
        return identifier in self._models

    def add_model(self, identifier: str, model: M):
        self._models[identifier] = model
        self._link_model(model)

    def add_models(self, models: dict[str, M]):
        for identifier, model in models.items():
            self.add_model(identifier, model)

    def get_model(self, identifier: str) -> M:
        keys = identifier.split(".", 1)
        if self.has_model(keys[0]):
            if len(keys) == 1:
                return self._models[identifier]
            else:
                sub_model = self._models[keys[0]]
                if isinstance(sub_model, HasModels):
                    return sub_model.get_model(keys[1])
                raise TypeError(
                    f"Model with identifier '{identifier}' does not have sub-models."
                )
        raise ValueError(f"Model with identifier '{identifier}' not found.")

    def get_models(self) -> t.Iterable[tuple[str, M]]:
        return self._models.items()

    def _link_model(self, model: M):
        tl.dlink(
            (self, "locked"),
            (model, "locked"),
        )
        tl.dlink(
            (model, "blockers"),
            (self, "blockers"),
        )
        for dependency in model.dependencies:
            dependency_parts = dependency.rsplit(".", 1)
            if len(dependency_parts) == 1:  # from parent
                target_model = self
                trait = dependency
            else:  # from sibling
                sibling, trait = dependency_parts
                target_model = self.get_model(sibling)
            tl.dlink(
                (target_model, trait),
                (model, trait),
            )


class HasProcess(tl.HasTraits):
    process_uuid = tl.Unicode(None, allow_none=True)
    monitor_counter = tl.Int(0)  # used for continuous updates

    @property
    def process(self) -> orm.WorkChainNode | None:
        if not self.process_uuid:
            return None
        try:
            return t.cast(orm.WorkChainNode, orm.load_node(self.process_uuid))
        except NotExistent:
            return None

    @property
    def has_process(self):
        return self.process is not None

    @property
    def inputs(self) -> orm.NodeLinksManager | None:
        return self.process.inputs if self.has_process else None

    @property
    def properties(self) -> list:
        return (
            self.inputs.properties.base.attributes.get("list")
            if self.inputs and "properties" in self.inputs
            else []
        )

    @property
    def outputs(self) -> orm.NodeLinksManager | None:
        return self.process.outputs if self.has_process else None


class Confirmable(tl.HasTraits):
    confirmed = tl.Bool(False)

    confirmation_exceptions = [
        "confirmed",
    ]

    def confirm(self):
        self.confirmed = True

    @tl.observe(tl.All)
    def _on_any_change(self, change):
        if change and change["name"] not in self.confirmation_exceptions:
            self._unconfirm()

    def _unconfirm(self):
        self.confirmed = False


class HasBlockers(tl.HasTraits):
    blockers = tl.List(tl.Unicode())

    @property
    def is_blocked(self):
        return any(self.blockers)

    def update_blockers(self):
        blockers = list(self._check_blockers() or [])
        if isinstance(self, HasModels):
            for _, model in self.get_models():
                if isinstance(model, HasBlockers):
                    blockers.extend(model.blockers)
        self.blockers = blockers

    def reset_blockers(self):
        self.blockers = []
        if isinstance(self, HasModels):
            for _, model in self.get_models():
                if isinstance(model, HasBlockers):
                    model.reset_blockers()

    def _check_blockers(self) -> t.Generator[str, None, None] | None:
        pass
