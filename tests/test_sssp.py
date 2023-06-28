def test_setup_codes(sssp):
    from aiidalab_qe.app.sssp import SSSPInstallWidget, pseudos_to_install

    widget = SSSPInstallWidget(auto_start=False)
    widget.set_message("abc")

    assert len(pseudos_to_install()) == 0
