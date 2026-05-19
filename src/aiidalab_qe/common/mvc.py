from __future__ import annotations

import typing as t

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

    Attributes
    ----------
    `identifier` : `str`
        The identifier for this model.
    `dependencies` : `list[str]`
        A list of dependencies for this model.
    `updated` : `bool`
        Whether the model has been updated.
    `locked` : `bool`
        Whether the model is locked.
    """

    identifier = ""
    dependencies: list[str] = []
    updated = False

    locked = tl.Bool(False)

    def __init__(self, state: dict | None = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._defaults = {}
        if state:
            self.set_model_state(state)

    def lock(self):
        """Lock the model, preventing any further changes."""
        self.locked = True

    def update(self, specific=""):
        """Updates the model.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of update.
        """
        pass

    def get_model_state(self) -> dict:
        """Retrieve the current model state."""
        raise NotImplementedError()

    def set_model_state(self, state: dict):
        """Preload the model state."""
        raise NotImplementedError()

    def reset(self):
        """Reset the model to present defaults."""
        pass

    def _get_default(self, trait: dict) -> t.Any:
        return self._defaults.get(trait, self.trait_defaults()[trait])
