"""Persist plugin activations that still require an AiiDA daemon restart."""

from pathlib import Path

from packaging.utils import canonicalize_name

from aiida.manage.configuration import get_profile
from aiidalab_qe.config import CONFIG_DIR_PATH


def _pending_restarts_dir() -> Path:
    return CONFIG_DIR_PATH / "pending-plugin-restarts"


def _restart_marker(plugin_name: str) -> Path:
    profile = get_profile()
    profile_id = profile.uuid if profile else "default"
    return _pending_restarts_dir() / profile_id / canonicalize_name(plugin_name)


def get_restart_failure(plugin_name: str) -> str | None:
    """Return an activation warning while this profile needs a daemon restart."""
    if _restart_marker(plugin_name).exists():
        return "Plugin activation requires a successful AiiDA daemon restart."
    return None


def set_restart_pending(plugin_name: str) -> None:
    """Persist the need to restart before attempting the daemon restart."""
    marker = _restart_marker(plugin_name)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.touch()


def clear_restart_pending(plugin_name: str) -> None:
    """Clear pending activations for this profile after a successful restart."""
    for marker in _restart_marker(plugin_name).parent.glob("*"):
        marker.unlink(missing_ok=True)
