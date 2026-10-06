def test_entry_point():
    from aiidalab_qe.plugins.utils import get_entries

    entries = get_entries()
    assert "bands" in entries
    assert "pdos" in entries
    assert "workchain" in entries["bands"]

    entries_list = list(entries.keys())
    prioritized_entries = ["electronic_structure", "bands", "pdos"]
    for i, prioritized_entry in enumerate(prioritized_entries):
        assert entries_list.index(prioritized_entry) == i, (
            f"Entry point {prioritized_entry} is not in the expected position."
        )


def test_entry_point_filter_skips_loading(monkeypatch):
    import importlib_metadata

    from aiidalab_qe.plugins.utils import get_entries

    class EntryPoint:
        name = "incompatible"

        def load(self):
            raise AssertionError("filtered entry point must not be loaded")

    monkeypatch.setattr(importlib_metadata, "entry_points", lambda **_: [EntryPoint()])

    entries = get_entries(
        "test.group",
        entry_point_filter=lambda _entry_point: False,
    )

    assert entries == {}


def test_failed_activation_filters_entry_point_unless_explicitly_included(
    monkeypatch, tmp_path
):
    import importlib_metadata

    from aiidalab_qe.plugins import state
    from aiidalab_qe.plugins.utils import get_entries

    class Distribution:
        name = "my-plugin"

    class EntryPoint:
        name = "plugin-entry"
        dist = Distribution()

        def load(self):
            return {"loaded": True}

    monkeypatch.setattr(
        state, "ACTIVATION_STATE_PATH", tmp_path / "plugin-activation.json"
    )
    state.set_activation_failure("my_plugin", "validation failed")
    monkeypatch.setattr(importlib_metadata, "entry_points", lambda **_: [EntryPoint()])

    assert get_entries("test.group") == {}
    assert get_entries("test.group", include_activation_failures=True) == {
        "plugin-entry": {"loaded": True}
    }
