def test_setup_codes(pw_code, projwfc_code, dos_code):
    from aiidalab_qe.app.setup_codes import QESetupWidget, codes_are_setup

    widget = QESetupWidget(auto_start=False)
    widget.set_message()

    assert codes_are_setup() is True
