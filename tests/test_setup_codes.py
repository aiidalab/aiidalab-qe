from aiidalab_qe.setup import codes


def test_code_setup_check_loads_profile(monkeypatch):
    calls = []

    def load_profile():
        calls.append(None)

    def load_code(identifier):
        assert identifier == f"pw-{codes.QE_VERSION}@localhost"

    monkeypatch.setattr(codes, "CODE_NAMES", ("pw",))
    monkeypatch.setattr(codes, "load_profile", load_profile)
    monkeypatch.setattr(codes, "load_code", load_code)

    assert codes.codes_are_setup("localhost")
    assert calls == [None]
