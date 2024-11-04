import typing as t

import traitlets as tl

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

T = t.TypeVar("T")


class MissingMixinError(Exception):
    """Raised when no mixin is found in a class definition."""

    def __init__(self, *args: object) -> None:
        super().__init__("expected at least one mixin in class definion", *args)


class MetaMixinTraitsRegister(tl.MetaHasTraits):
    """A metaclass to register traits from `HasTraits`-subclassed mixins.

    The metaclass removes the `HasTraits` base class from the mixin to
    avoid MRO conflicts.
    """

    def __new__(cls, name, bases, classdict):
        if len(bases) == 1 and not issubclass(bases[0], tl.HasTraits):
            raise MissingMixinError()
        for base in bases:
            if issubclass(base, tl.HasTraits):
                for name, trait in base.class_traits().items():
                    if name not in classdict:
                        classdict[name] = trait
        bases = tuple(filter(lambda base: base is not tl.HasTraits, bases))
        return super().__new__(cls, name, bases, classdict)


class HasTraitsAndMixins(tl.HasTraits, metaclass=MetaMixinTraitsRegister):
    """An extension of `traitlet`'s `HasTraits` to support trait-ful mixins."""


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


class Confirmable(tl.HasTraits):
    confirmed = tl.Bool(False)

    def confirm(self):
        self.confirmed = True

    @tl.observe(tl.All)
    def unconfirm(self, change=None):
        if change and change["name"] != "confirmed":
            self.confirmed = False
