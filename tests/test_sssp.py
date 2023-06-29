def test_sssp(sssp):
    """Test the setup of SSSP pseudopotentials."""
    from aiidalab_qe.app.sssp import install, install_pseudos, pseudos_to_install

    assert len(pseudos_to_install()) == 0
    #
    install()
    #
    install_pseudos(["SSSP/1.2/PBE/efficiency"])


def test_sssp_widget(sssp):
    """Test the widget for installing SSSP pseudopotentials.""" ""
    from aiidalab_qe.app.sssp import SSSPInstallWidget

    widget = SSSPInstallWidget(auto_start=False)
    widget.set_message("abc")
