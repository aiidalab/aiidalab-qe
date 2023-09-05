def test_entry_point():
    from aiidalab_qe.app.utils import get_entries

    entries = get_entries()
    print(entries)
    assert "bands" in entries
