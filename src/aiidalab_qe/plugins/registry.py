"""Plugin registry data and live installation status."""

from __future__ import annotations

import asyncio
import logging
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import metadata
from importlib.resources import files

import yaml
from packaging.requirements import InvalidRequirement, Requirement
from packaging.specifiers import InvalidSpecifier, SpecifierSet
from packaging.utils import canonicalize_name
from packaging.version import InvalidVersion, Version

from aiidalab_qe import __version__

LOGGER = logging.getLogger(__name__)
POST_INSTALL_VERIFY_TIMEOUT = 60


def poll_post_install(
    plugin_id: str,
    post_install: str | None,
    profile: str = "default",
) -> tuple[bool | None, str | None]:
    """Check post-install state using the plugin's read-only verifier."""
    if not post_install:
        return None, None

    command = [
        sys.executable,
        "-m",
        plugin_id.replace("-", "_"),
        "verify-post-install",
        "--profile",
        profile,
    ]
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=POST_INSTALL_VERIFY_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return False, (
            "Post-install verification timed out after "
            f"{POST_INSTALL_VERIFY_TIMEOUT} seconds."
        )
    except OSError as error:
        LOGGER.warning(
            "Could not start post-install verification for %s: %s", plugin_id, error
        )
        return False, f"Could not start post-install verification: {error}"

    if result.returncode == 0:
        return True, None

    details = (result.stderr or result.stdout or "").strip()
    if len(details) > 500:
        details = f"{details[:497]}..."
    message = "Post-install verification failed."
    if details:
        message = f"{message} {details}"
    return False, message


def poll_plugin_registration(plugin_id: str) -> str | None:
    """Check entry-point loading and app integration for an installed plugin."""
    from aiidalab_qe.plugins.utils import test_plugin_functionality

    registered, error = test_plugin_functionality(plugin_id)
    if registered:
        return None
    if not registered and not error:
        error = "Plugin entry-point validation failed without an error message."
    return error


def poll_daemon_status() -> str | None:
    """Return a warning when the current AiiDA daemon is unavailable."""
    from aiida.engine.daemon.client import get_daemon_client

    try:
        if get_daemon_client().is_daemon_running:
            return None
    except Exception as error:
        LOGGER.warning("Could not determine AiiDA daemon status: %s", error)
        return f"Could not determine AiiDA daemon status: {error}"

    return "The AiiDA daemon is not running. Restart it, then refresh plugin status."


@dataclass(frozen=True)
class QeAppPluginData:
    title: str
    description: str
    pip: str | None = None
    github: str | None = None
    author: str = "N/A"
    documentation: str | None = None
    post_install: str | None = None
    status: str = ""
    category: str = "calculation"
    requires_aiidalab_qe: str | None = None
    extra: dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_mapping(cls, plugin_name: str, data: Mapping) -> QeAppPluginData:
        """Create a plugin record after validating fields consumed by the app."""
        if not isinstance(plugin_name, str) or not plugin_name.strip():
            raise ValueError("plugin registry key must be a non-empty string")
        if not isinstance(data, Mapping):
            raise TypeError(f"{plugin_name}: plugin entry must be a mapping")

        def required_string(field: str) -> str:
            value = data.get(field)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{plugin_name}: '{field}' must be a non-empty string")
            return value

        def optional_string(field: str, default: str | None = None) -> str | None:
            value = data.get(field, default)
            if value is None:
                return None
            if not isinstance(value, str):
                raise TypeError(f"{plugin_name}: '{field}' must be a string")
            return value or None

        title = required_string("title")
        description = required_string("description")

        pip_requirement = optional_string("pip")
        github = optional_string("github")
        if not pip_requirement and not github:
            raise ValueError(f"{plugin_name}: either 'pip' or 'github' is required")

        try:
            Requirement(pip_requirement or plugin_name)
        except InvalidRequirement as error:
            field_name = "pip" if pip_requirement else "top-level key"
            raise ValueError(
                f"{plugin_name}: invalid '{field_name}' requirement: {error}"
            ) from error

        app_requirement = optional_string("requires_aiidalab_qe")
        if app_requirement:
            try:
                SpecifierSet(app_requirement)
            except InvalidSpecifier as error:
                raise ValueError(
                    f"{plugin_name}: invalid 'requires_aiidalab_qe' specifier: {error}"
                ) from error

        known_fields = {
            "title",
            "description",
            "pip",
            "github",
            "author",
            "documentation",
            "post_install",
            "status",
            "category",
            "requires_aiidalab_qe",
        }

        return cls(
            title=title,
            description=description,
            pip=pip_requirement,
            github=github,
            author=optional_string("author", "N/A") or "N/A",
            documentation=optional_string("documentation"),
            post_install=optional_string("post_install"),
            status=optional_string("status", "") or "",
            category=optional_string("category", "calculation") or "calculation",
            requires_aiidalab_qe=app_requirement,
            extra={
                key: value for key, value in data.items() if key not in known_fields
            },
        )

    def as_mapping(self) -> dict:
        """Expose normalized metadata for existing registry consumers."""
        return {
            **self.extra,
            "title": self.title,
            "description": self.description,
            "pip": self.pip,
            "github": self.github,
            "author": self.author,
            "documentation": self.documentation,
            "post_install": self.post_install,
            "status": self.status,
            "category": self.category,
            "requires_aiidalab_qe": self.requires_aiidalab_qe,
        }


def load_plugin_config() -> Mapping:
    """Load raw plugin definitions from the packaged registry."""
    try:
        content = (
            files("aiidalab_qe.plugins")
            .joinpath("plugins.yaml")
            .read_text(encoding="utf-8")
        )
        data = yaml.safe_load(content)
    except (OSError, yaml.YAMLError) as error:
        raise ValueError(f"Could not read packaged plugin registry: {error}") from error

    if data is None:
        return {}
    if not isinstance(data, Mapping):
        raise TypeError("plugin registry must be a mapping of names to entries")
    return data


class PluginRegistry:
    """Validated, data-only view of the packaged plugin registry."""

    def __init__(
        self,
        plugins: Mapping[str, QeAppPluginData] | None = None,
        config_error: str | None = None,
    ):
        self.plugins = dict(plugins or {})
        self.config_error = config_error

    @classmethod
    def load(cls) -> PluginRegistry:
        """Load and validate a registry, retaining errors for UI presentation."""
        try:
            raw_plugins = load_plugin_config()
            plugins = {
                name: QeAppPluginData.from_mapping(name, data)
                for name, data in raw_plugins.items()
            }
        except (
            OSError,
            TypeError,
            ValueError,
            yaml.YAMLError,
        ) as error:
            LOGGER.exception("Could not load packaged plugin registry")
            return cls(config_error=str(error))

        return cls(plugins=plugins)

    @property
    def data(self) -> dict:
        """Return normalized registry entries for existing consumers."""
        return {
            name: plugin_data.as_mapping() for name, plugin_data in self.plugins.items()
        }


@lru_cache(maxsize=1)
def get_plugin_registry() -> PluginRegistry:
    """Return the packaged registry cached for this Python process."""
    return PluginRegistry.load()


def is_version_compatible(required_version: str | None) -> bool:
    """Check whether the installed app version satisfies a requirement."""
    if not required_version:
        return True

    try:
        specifier = SpecifierSet(required_version)
        installed_version = Version(__version__)
        return specifier.contains(installed_version, prereleases=True)
    except (InvalidSpecifier, InvalidVersion) as error:
        LOGGER.warning(
            "Could not parse version requirement %r: %s",
            required_version,
            error,
        )
        return False


def is_plugin_installed(plugin_name: str) -> bool:
    """Check whether a package is installed."""
    try:
        metadata.version(plugin_name)
    except metadata.PackageNotFoundError:
        return False
    return True


def get_plugin_version_info(
    plugin_name: str,
    pip_requirement: str | None,
) -> tuple[str | None, bool, str | None]:
    """Return the installed version, compatibility, and any requirement error."""
    try:
        requirement = Requirement(pip_requirement or plugin_name)
    except InvalidRequirement as error:
        return None, False, str(error)

    try:
        installed_version = metadata.version(requirement.name)
    except metadata.PackageNotFoundError:
        return None, True, None

    try:
        version = Version(installed_version)
    except InvalidVersion as error:
        return installed_version, False, str(error)

    compatible = not requirement.specifier or requirement.specifier.contains(
        version,
        prereleases=True,
    )
    return installed_version, compatible, None


@dataclass(frozen=True)
class PluginStatus:
    installed_version: str | None
    plugin_compatible: bool
    requirement_error: str | None
    app_compatible: bool
    registration_error: str | None
    post_install_verified: bool | None
    post_install_error: str | None
    is_installed: bool
    distribution_name: str
    requirement: str

    @property
    def is_incompatible(self) -> bool:
        return self.is_installed and bool(
            not self.plugin_compatible
            or self.requirement_error
            or not self.app_compatible
            or self.registration_error
            or self.post_install_verified is False
        )


def get_plugin_status(
    plugin_name: str,
    plugin_data: QeAppPluginData,
    *,
    check_runtime: bool = True,
) -> PluginStatus:
    """Collect live compatibility and activation status for one plugin."""
    requirement_text = plugin_data.pip or plugin_name
    registration_error = None
    post_install_verified = None
    post_install_error = None

    (
        installed_version,
        plugin_compatible,
        requirement_error,
    ) = get_plugin_version_info(plugin_name, requirement_text)

    app_compatible = is_version_compatible(plugin_data.requires_aiidalab_qe)
    is_installed = installed_version is not None

    try:
        requirement = Requirement(requirement_text)
    except InvalidRequirement:
        distribution_name = canonicalize_name(plugin_name)
        is_installed = is_plugin_installed(plugin_name)
    else:
        distribution_name = canonicalize_name(requirement.name)
        requirement_text = str(requirement)

    post_install = getattr(plugin_data, "post_install", None)
    if (
        is_installed
        and plugin_compatible
        and app_compatible
        and not requirement_error
        and check_runtime
    ):

        def verify_runtime():
            registration = poll_plugin_registration(plugin_name)
            if post_install:
                verified, error = poll_post_install(plugin_name, post_install)
            else:
                verified, error = None, None
            return registration, verified, error

        (
            registration_error,
            post_install_verified,
            post_install_error,
        ) = verify_runtime()

    return PluginStatus(
        installed_version=installed_version,
        plugin_compatible=plugin_compatible,
        requirement_error=requirement_error,
        app_compatible=app_compatible,
        registration_error=registration_error,
        post_install_verified=post_install_verified,
        post_install_error=post_install_error,
        is_installed=is_installed,
        distribution_name=distribution_name,
        requirement=requirement_text,
    )


async def get_plugin_statuses(
    plugins: Mapping[str, QeAppPluginData],
    *,
    check_runtime: bool = True,
) -> dict[str, PluginStatus]:
    """Collect multiple plugin statuses concurrently without blocking the UI loop."""
    names = list(plugins)
    statuses = await asyncio.gather(
        *(
            asyncio.to_thread(
                get_plugin_status,
                name,
                plugins[name],
                check_runtime=check_runtime,
            )
            for name in names
        )
    )
    return dict(zip(names, statuses, strict=True))
