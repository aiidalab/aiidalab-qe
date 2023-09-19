def test_entry_point():
    from aiidalab_qe.app.utils import get_entries

    entries = get_entries()
    assert "bands" in entries
    assert "pdos" in entries
    assert "workchain" in entries["bands"]
