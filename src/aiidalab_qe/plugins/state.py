"""Persist plugin activation failures across app and kernel reloads."""

import json
import logging
from pathlib import Path

from packaging.utils import canonicalize_name

LOGGER = logging.getLogger(__name__)

ACTIVATION_STATE_PATH = (
    Path.home() / ".aiidalab" / "quantum-espresso" / "plugin-activation.json"
)


def _read_state() -> dict[str, object] | None:
    """Read the activation state, returning None when it is invalid or unreadable."""
    try:
        state = json.loads(ACTIVATION_STATE_PATH.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    except (OSError, json.JSONDecodeError) as error:
        LOGGER.warning("Could not read plugin activation state: %s", error)
        return None

    if not isinstance(state, dict):
        LOGGER.warning("Plugin activation state must be a JSON object")
        return None

    return state


def _write_state(state: dict[str, object]) -> None:
    """Write activation state using an atomic file replacement."""
    ACTIVATION_STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = ACTIVATION_STATE_PATH.with_suffix(".tmp")
    temporary_path.write_text(json.dumps(state, indent=2), encoding="utf-8")
    temporary_path.replace(ACTIVATION_STATE_PATH)


def get_activation_failure(plugin_id: str) -> str | None:
    """Return the saved activation failure for a plugin, if present."""
    plugin_id = canonicalize_name(plugin_id)
    state = _read_state()

    if state is None:
        return None

    failure = state.get(plugin_id)
    return failure if isinstance(failure, str) else None


def set_activation_failure(plugin_id: str, message: str) -> None:
    """Persist a plugin activation failure using an atomic file replacement."""
    plugin_id = canonicalize_name(plugin_id)
    state = _read_state() or {}
    state[plugin_id] = message
    _write_state(state)


def clear_activation_failure(plugin_id: str) -> None:
    """Clear a plugin activation failure, preserving other plugins' state."""
    plugin_id = canonicalize_name(plugin_id)
    state = _read_state()

    if not state or plugin_id not in state:
        return

    del state[plugin_id]

    if state:
        _write_state(state)
    else:
        ACTIVATION_STATE_PATH.unlink(missing_ok=True)
