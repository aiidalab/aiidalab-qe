import typing as t

import traitlets as tl

from aiida import orm
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
    process_node = tl.Instance(orm.Node, allow_none=True)

    @property
    def has_process(self):
        return self.process_node is not None


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
