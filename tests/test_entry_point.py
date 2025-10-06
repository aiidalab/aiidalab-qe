def test_entry_point():
    from aiidalab_qe.app.utils import get_entries

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
