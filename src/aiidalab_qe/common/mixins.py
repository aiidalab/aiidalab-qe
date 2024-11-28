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


class HasModels(t.Generic[T]):
    def __init__(self):
        self._models: dict[str, T] = {}

    def has_model(self, identifier):
        return identifier in self._models

    def add_model(self, identifier, model):
        self._models[identifier] = model
        self._link_model(model)

    def add_models(self, models: dict[str, T]):
        for identifier, model in models.items():
            self.add_model(identifier, model)

    def get_model(self, identifier) -> T:
        if self.has_model(identifier):
            return self._models[identifier]
        raise ValueError(f"Model with identifier '{identifier}' not found.")

    def get_models(self) -> t.Iterable[tuple[str, T]]:
        return self._models.items()

    def _link_model(self, model: T):
        raise NotImplementedError()


class HasProcess(tl.HasTraits):
    process_uuid = tl.Unicode(allow_none=True)
    monitor_counter = tl.Int(0)  # used for continuous updates

    @property
    def inputs(self):
        process_node = self.fetch_process_node()
        return process_node.inputs if process_node else []

    @property
    def properties(self):
        process_node = self.fetch_process_node()
        return process_node.inputs.properties if process_node else []

    @property
    def outputs(self):
        process_node = self.fetch_process_node()
        return process_node.outputs if process_node else []

    def fetch_process_node(self):
        try:
            return orm.load_node(self.process_uuid) if self.process_uuid else None
        except NotExistent:
            return None


class Confirmable(tl.HasTraits):
    confirmed = tl.Bool(False)

    def confirm(self):
        self.confirmed = True

    @tl.observe(tl.All)
    def _on_any_change(self, change):
        if change and change["name"] != "confirmed":
            self._unconfirm()

    def _unconfirm(self):
        self.confirmed = False
