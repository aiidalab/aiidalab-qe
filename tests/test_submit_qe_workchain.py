def test_code_not_selected(submit_app_generator):
    """Test if there is an error when the code is not selected."""
    app = submit_app_generator()
    app.submit_step.dos_code.value = None
    app.submit_step._create_builder()


def test_reload_selected_code(submit_app_generator):
    """Test set_selected_codes method."""
    from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep

    app = submit_app_generator()
    submit_step = app.submit_step

    submit_step._create_builder()

    new_submit_step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    new_submit_step.set_selected_codes(submit_step.ui_parameters["codes"])

    assert new_submit_step.pw_code.value == submit_step.pw_code.value
    assert new_submit_step.dos_code.value == submit_step.dos_code.value
    assert new_submit_step.projwfc_code.value == submit_step.projwfc_code.value


def test_create_builder_default(
    data_regression,
    submit_app_generator,
):
    """ "Test the creation of the workchain builder.

    metal, non-magnetic
    """

    app = submit_app_generator(properties=["bands", "pdos"])
    submit_step = app.submit_step

    builder = submit_step._create_builder()

    # check and validate the builder
    got = builder_to_readable_dict(builder)

    # regression test
    data_regression.check(got)


def test_create_builder_insulator(
    submit_app_generator,
):
    """ "Test the creation of the workchain builder.

    insulator, non-magnetic, no smearing
    the occupation type is set to fixed, smearing and degauss should not be set"""

    app = submit_app_generator(
        electronic_type="insulator", properties=["bands", "pdos"]
    )
    submit_step = app.submit_step

    builder = submit_step._create_builder()

    # check and validate the builder
    got = builder_to_readable_dict(builder)

    assert got["bands"]["scf"]["pw"]["parameters"]["SYSTEM"]["occupations"] == "fixed"
    assert "smearing" not in got["bands"]["scf"]["pw"]["parameters"]["SYSTEM"]


def test_create_builder_advanced_settings(
    submit_app_generator,
):
    """Test the creation of the workchain builder with advanced settings

    -metal
    -collinear
    -tot_charge
    -initial_magnetic_moments
    """

    app = submit_app_generator(
        electronic_type="metal",
        spin_type="collinear",
        tot_charge=1.0,
        initial_magnetic_moments=0.1,
        properties=["bands", "pdos"],
    )
    submit_step = app.submit_step

    builder = submit_step._create_builder()

    # check and validate the builder
    got = builder_to_readable_dict(builder)

    # test tot_charge is updated in the three steps
    for parameters in [
        got["relax"]["base"],
        got["bands"]["scf"],
        got["pdos"]["scf"],
        got["pdos"]["nscf"],
    ]:
        assert parameters["pw"]["parameters"]["SYSTEM"]["tot_charge"] == 1.0

    # test initial_magnetic_moments set 'starting_magnetization' in pw.in
    assert (
        got["relax"]["base"]["pw"]["parameters"]["SYSTEM"]["starting_magnetization"][
            "Si"
        ]
        == 0.025
    )


def builder_to_readable_dict(builder):
    """transverse the builder and return a dictionary with readable values."""
    from aiida import orm
    from aiida.engine import ProcessBuilderNamespace
    from aiida.plugins import DataFactory

    UpfData = DataFactory("pseudo.upf")

    ignore_keys = ["metadata", "monitors", "code", "structure"]

    readable_dict = {}
    for k, v in builder.items():
        if k in ignore_keys:
            continue
        if isinstance(v, UpfData):
            readable_dict[k] = v.filename
        elif isinstance(v, (dict, ProcessBuilderNamespace)):
            readable_dict[k] = builder_to_readable_dict(v)
        elif isinstance(v, orm.Dict):
            readable_dict[k] = v.get_dict()
        elif isinstance(v, (orm.Int, orm.Float, orm.Str, orm.Bool)):
            readable_dict[k] = v.value
        elif isinstance(v, orm.List):
            readable_dict[k] = v.get_list()
        else:
            readable_dict[k] = v

    return readable_dict
