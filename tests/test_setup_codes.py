def test_setup_codes(pw_code, projwfc_code, dos_code, fixture_localhost):
    from aiida.orm import load_code
    from aiida.tools import delete_nodes

    from aiidalab_qe.app.setup_codes import (
        QE_VERSION,
        _code_is_setup,
        _setup_code,
        codes_are_setup,
    )

    assert codes_are_setup() is True
    assert _code_is_setup("pw") is True
    # assert qe_installed() is False
    # delete codes if they exist
    try:
        pw_code = load_code(f"pw-{QE_VERSION}@localhost")
        print("pw_code: ", pw_code)
        delete_nodes([pw_code.pk], dry_run=False)
    except Exception:
        pass
    assert _code_is_setup("pw") is False
    _setup_code("pw")
    assert _code_is_setup("pw") is True


def test_widget(pw_code, projwfc_code, dos_code):
    from aiidalab_qe.app.setup_codes import QESetupWidget

    widget = QESetupWidget(auto_start=False)
    widget.set_message("abc")
