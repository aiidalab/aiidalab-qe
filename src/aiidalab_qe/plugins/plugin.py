"""Plugin interface exposed to third-party AiiDAlab QE plugins."""

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Any


class QeAppPlugin(ABC):
    """Base class for a plugin's lazily loaded app integrations.

    One instance is cached for the app process. Keep construction lightweight and
    avoid storing per-calculation state here. Put feature imports inside methods so
    features load only when the app requests them.
    """

    @property
    @abstractmethod
    def outline(self) -> str:
        """Return the title used for the plugin's Step 2 outline row."""

    @property
    def description(self) -> str:
        """Return optional outline help text."""
        return ""

    def structure_examples(self) -> Any | None:
        """Return structure examples, if available."""
        return None

    def importers(self) -> tuple[Callable[[], Any], ...]:
        """Return importer component classes contributed by this plugin."""
        return ()

    def editors(self) -> tuple[Callable[[], Any], ...]:
        """Return editor component classes contributed by this plugin."""
        return ()

    def configuration(self) -> Any | None:
        """Return configuration panel and model, if the plugin has this feature."""
        return None

    def resources(self) -> Any | None:
        """Return resource panel and model, if the plugin has this feature."""
        return None

    def workchain(self) -> Any | None:
        """Return workchain configuration, if the plugin has this feature."""
        return None

    def result(self) -> Any | None:
        """Return a mapping of the plugin's result step component and model."""
        return None

    def guides(self) -> Any | None:
        """Return a mapping of the plugin's guide items."""
        return None

    def metadata(self) -> Any | None:
        """Return metadata for the plugin."""
        return None

    def examples(self) -> Any | None:
        """Return a custom example-archive source configuration, if available."""
        return None

    def load_feature(self, name: str) -> Any | None:
        """Load a named optional plugin feature."""
        if name == "outline":
            return self.outline
        if name == "description":
            return self.description
        if name not in {
            "structure_examples",
            "importers",
            "editors",
            "configuration",
            "resources",
            "workchain",
            "result",
            "guides",
            "metadata",
            "examples",
        }:
            raise KeyError(name)
        return getattr(self, name)()
