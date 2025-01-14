import traitlets as tl


class MetaHasTraitsLast(tl.MetaHasTraits):
    """A metaclass that ensures that `HasTraits` is pushed to the end of the MRO.

    This is useful when injecting trait-ful mixins into the MRO of a class."""

    def __new__(cls, name, bases, classdict):
        bases = (
            *filter(lambda base: base is not tl.HasTraits, bases),
            tl.HasTraits,
        )
        return super().__new__(cls, name, bases, classdict)


class Model(tl.HasTraits, metaclass=MetaHasTraitsLast):
    """A parent class for all MVC models.

    The class extends `traitlet`'s `HasTraits` and uses a metaclass to
    ensure that the MRO is reordered such that `HasTraits`, which may
    be extended by other bases (mixins), is pushed to the end of the
    MRO to avoid conflicts - see C3 linearization.
    """
