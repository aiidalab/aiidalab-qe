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
