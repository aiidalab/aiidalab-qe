"""Discover cached plugin instances and load their optional features."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from functools import lru_cache
from importlib import invalidate_caches
from types import MappingProxyType

from importlib_metadata import distributions, entry_points

from aiidalab_qe.plugins import QeAppPlugin

PLUGIN_ENTRY_POINT_GROUP = "aiidalab_qe.plugins"
PRIORITY_PLUGINS = ("electronic_structure", "bands", "pdos")


def print_error(plugin, error, feature_name: str | None = None) -> None:
    plugin_name = getattr(plugin, "name", plugin)
    if feature_name:
        message = f"plugin feature {feature_name} for {plugin_name}"
    else:
        message = f"plugin {plugin_name}"
    print(f"\033[91mFailed to load {message}.\033[0m")
    print(
        "\033[93mThis may be due to compatibility issues with the current QEApp version.\033[0m"
    )
    print("\033[93mPlease contact the plugin author for further assistance.\033[0m")
    print(
        "\033[93mThus, the plugin will not be available. However, you can still use the rest of the app.\033[0m"
    )
    print(f"\033[91mError message: {error}\033[0m\n")


@dataclass(frozen=True)
class PluginCatalog:
    """Plugin instances discovered once for the current Python app process."""

    plugins: Mapping[str, QeAppPlugin]

    def get_features(
        self,
        feature_name: str,
        plugin_filter: Callable[[str, QeAppPlugin], bool] | None = None,
    ) -> dict[str, object]:
        features = {}
        for name, plugin in self.plugins.items():
            if plugin_filter is not None and not plugin_filter(name, plugin):
                continue
            try:
                feature = plugin.load_feature(feature_name)
            except Exception as error:
                print_error(name, error, feature_name=feature_name)
                continue
            if feature:
                features[name] = feature
        return features


@lru_cache(maxsize=1)
def load_plugin_catalog() -> PluginCatalog:
    """Discover plugin classes and cache one instance of each."""
    sorted_entry_points = sorted(
        entry_points(group=PLUGIN_ENTRY_POINT_GROUP),
        key=lambda entry_point: (
            entry_point.name not in PRIORITY_PLUGINS,
            PRIORITY_PLUGINS.index(entry_point.name)
            if entry_point.name in PRIORITY_PLUGINS
            else entry_point.name,
        ),
    )

    plugins = {}
    for entry_point in sorted_entry_points:
        if entry_point.name in plugins:
            continue
        try:
            plugin_class = entry_point.load()
            if not isinstance(plugin_class, type) or not issubclass(
                plugin_class,
                QeAppPlugin,
            ):
                print_error(
                    entry_point,
                    TypeError("entry point must reference a QeAppPlugin subclass"),
                )
                continue
            plugins[entry_point.name] = plugin_class()
        except Exception as error:
            print_error(entry_point, error)

    return PluginCatalog(MappingProxyType(plugins))


def get_plugins() -> Mapping[str, QeAppPlugin]:
    """Return all cached plugin instances in registration order."""
    return load_plugin_catalog().plugins


def get_plugin_feature(plugin: QeAppPlugin, feature_name: str):
    """Load one feature from a plugin instance."""
    return plugin.load_feature(feature_name)


def get_plugin_features(
    feature_name: str,
    plugin_filter: Callable[[str, QeAppPlugin], bool] | None = None,
) -> dict[str, object]:
    """Load one feature from each selected plugin."""
    return load_plugin_catalog().get_features(feature_name, plugin_filter)


def get_plugin_entry_points_for_package(
    package_name: str,
    group: str = PLUGIN_ENTRY_POINT_GROUP,
):
    """Find plugin entry points for an installed package."""
    distribution = next(
        (
            distribution
            for distribution in distributions()
            if distribution.metadata.get("Name") == package_name
        ),
        None,
    )
    if distribution is None:
        raise ValueError(f"Package '{package_name}' not found.")
    return [
        entry_point
        for entry_point in distribution.entry_points
        if entry_point.group == group
    ]


def test_plugin_functionality(plugin_name):
    """Validate the functionality of a QE app plugin."""
    try:
        invalidate_caches()
        for entry_point in get_plugin_entry_points_for_package(plugin_name):
            plugin_class = entry_point.load()
            if not isinstance(plugin_class, type) or not issubclass(
                plugin_class,
                QeAppPlugin,
            ):
                return (
                    False,
                    "Plugin entry point must reference a QeAppPlugin subclass.",
                )
    except Exception as error:
        return False, f"Failed to validate plugin {plugin_name}: {error}"
    return True, ""
