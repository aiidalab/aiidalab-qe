def test_codes(pw_code, dos_code, projwfc_code, pw_code_70):
    """Test get and set codes."""
    from aiida.orm import load_code

    from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
    from aiidalab_qe.app.submit import SubmitQeAppWorkChainStep

    step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    # get the default codes
    assert (
        load_code(step.pw_code.value).full_label
        == DEFAULT_PARAMETERS["codes"]["pw_code"]
    )
    assert (
        load_code(step.dos_code.value).full_label
        == DEFAULT_PARAMETERS["codes"]["dos_code"]
    )
    assert (
        load_code(step.projwfc_code.value).full_label
        == DEFAULT_PARAMETERS["codes"]["projwfc_code"]
    )
    # test `get_selected_codes` method
    codes = step.get_selected_codes()
    codes["pw_code"] = pw_code.uuid
    # set the codes manually
    step.pw_code.value = pw_code_70.uuid
    assert load_code(step.pw_code.value).full_label == pw_code_70.full_label
    # set `set_selected_codes` method
    step.set_selected_codes(parameters=DEFAULT_PARAMETERS)
    assert step.pw_code.value == pw_code.uuid


def test_resources():
    """Test get and set resources."""
    from aiidalab_qe.app.submit import SubmitQeAppWorkChainStep

    step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    step.resources_config.num_nodes.value = 2
    step.resources_config.num_cpus.value = 4
    step.parallelization.npools.value = 2
    assert step.resources_config.num_nodes.value == 2
    assert step.resources_config.num_cpus.value == 4
    assert step.parallelization.npools.value == 2
    # test `get_resource` method
    resources = step.get_resource()
    assert resources["num_machines"] == 2
    assert resources["num_mpiprocs_per_machine"] == 4
    assert resources["npool"] == 2


def test_create_builder_default(
    data_regression,
    submit_step_widget_generator,
):
    """ "Test the creation of the workchain builder.

    metal, non-magnetic
    """
    from aiidalab_qe.app.report import _generate_report_html

    submit_step = submit_step_widget_generator()

    builder = submit_step._create_builder()
    extra_parameters = submit_step._create_extra_report_parameters()

    # check and validate the builder
    got = builder_to_readable_dict(builder)

    # regression test
    data_regression.check(got)

    # test the report can be properly generated from the builder without errors
    builder_parameters = submit_step._extract_report_parameters(
        builder, extra_parameters
    )
    report_html = _generate_report_html(builder_parameters)

    # None in report_html means that the report not properly generated
    assert "None" not in report_html


def test_create_builder_insulator(
    submit_step_widget_generator,
):
    """ "Test the creation of the workchain builder.

    insulator, non-magnetic, no smearing
    the occupation type is set to fixed, smearing and degauss should not be set"""
    from aiidalab_qe.app.report import _generate_report_html

    submit_step = submit_step_widget_generator(
        electronic_type="insulator",
    )

    builder = submit_step._create_builder()
    extra_parameters = submit_step._create_extra_report_parameters()

    # check and validate the builder
    got = builder_to_readable_dict(builder)

    assert got["bands"]["scf"]["pw"]["parameters"]["SYSTEM"]["occupations"] == "fixed"
    assert "smearing" not in got["bands"]["scf"]["pw"]["parameters"]["SYSTEM"]

    # test the report can be properly generated from the builder without errors
    builder_parameters = submit_step._extract_report_parameters(
        builder, extra_parameters
    )
    report_html = _generate_report_html(builder_parameters)

    # None in report_html means that the report not properly generated
    assert "None" not in report_html


def builder_to_readable_dict(builder):
    """transverse the builder and return a dictionary with readable values."""
    from aiida import orm
    from aiida.engine import ProcessBuilderNamespace

    ignore_keys = ["metadata", "monitors", "pseudos", "code", "structure"]

    readable_dict = {}
    for k, v in builder.items():
        if k in ignore_keys:
            continue
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
