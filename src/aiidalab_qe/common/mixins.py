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

    def add_model(self, identifier, model):
        self._models[identifier] = model
        self._link_model(model)

    def get_model(self, identifier) -> T:
        if identifier in self._models:
            return self._models[identifier]
        raise ValueError(f"Model with identifier '{identifier}' not found.")

    def get_models(self) -> t.Iterable[tuple[str, T]]:
        return self._models.items()

    def _link_model(self, model: T):
        raise NotImplementedError()


class HasProcess(tl.HasTraits):
    process_uuid = tl.Unicode(allow_none=True)
    monitor_counter = tl.Int(0)

    process_node = None

    @property
    def has_process(self):
        return self.process_node is not None

    @property
    def inputs(self):
        return self.process_node.inputs if self.has_process else []

    @property
    def properties(self):
        return self.process_node.inputs.properties if self.has_process else []

    @property
    def outputs(self):
        return self.process_node.outputs if self.has_process else []

    def fetch_process_node(self):
        try:
            return orm.load_node(self.process_uuid) if self.process_uuid else None
        except NotExistent:
            return None

    @tl.observe("process_uuid")
    def _on_process_uuid_change(self, _):
        self.process_node = self.fetch_process_node()


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
