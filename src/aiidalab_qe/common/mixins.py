from __future__ import annotations

import typing as t

import traitlets as tl

from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

T = t.TypeVar("T")


class HasInputStructure(tl.HasTraits):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )

    @property
    def has_structure(self):
        return self.input_structure is not None

    @property
    def has_pbc(self):
        return not self.has_structure or any(self.input_structure.pbc)

    @property
    def has_tags(self):
        return any(
            not kind_name.isalpha()
            for kind_name in self.input_structure.get_kind_names()
        )


class HasModels(t.Generic[T]):
    def __init__(self):
        self._models: dict[str, T] = {}

    def has_model(self, identifier):
        return identifier in self._models

    def add_model(self, identifier, model: T):
        self._models[identifier] = model
        self._link_model(model)

    def add_models(self, models: dict[str, T]):
        for identifier, model in models.items():
            self.add_model(identifier, model)

    def get_model(self, identifier) -> T:
        keys = identifier.split(".", 1)
        if self.has_model(keys[0]):
            if len(keys) == 1:
                return self._models[identifier]
            else:
                return self._models[keys[0]].get_model(keys[1])
        raise ValueError(f"Model with identifier '{identifier}' not found.")

    def get_models(self) -> t.Iterable[tuple[str, T]]:
        return self._models.items()

    def _link_model(self, model: T):
        if not hasattr(model, "dependencies"):
            return
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
    def has_process(self):
        return self.fetch_process_node() is not None

    @property
    def inputs(self):
        process_node = self.fetch_process_node()
        return process_node.inputs if process_node else []

    @property
    def properties(self):
        process_node = self.fetch_process_node()
        # read the attributes directly instead of using the `get_list` method
        # to avoid error in case of the orm.List object being converted to a orm.Data object
        return (
            process_node.inputs.properties.base.attributes.get("list")
            if process_node
            else []
        )

    @property
    def outputs(self):
        process_node = self.fetch_process_node()
        return process_node.outputs if process_node else []

    def fetch_process_node(self) -> orm.ProcessNode | None:
        try:
            return orm.load_node(self.process_uuid) if self.process_uuid else None  # type: ignore
        except NotExistent:
            return None


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
    blocker_messages = tl.Unicode("")

    @property
    def is_blocked(self):
        return any(self.blockers)

    def update_blockers(self):
        blockers = list(self._check_blockers())
        if isinstance(self, HasModels):
            for _, model in self.get_models():
                if isinstance(model, HasBlockers):
                    blockers += model.blockers
        self.blockers = blockers

    def update_blocker_messages(self):
        if self.is_blocked:
            formatted = "\n".join(f"<li>{item}</li>" for item in self.blockers)
            self.blocker_messages = f"""
                <div class="alert alert-danger">
                    <b>The step is blocked due to the following reason(s):</b>
                    <ul>
                        {formatted}
                    </ul>
                </div>
            """
        else:
            self.blocker_messages = ""

    def _check_blockers(self):
        raise NotImplementedError
