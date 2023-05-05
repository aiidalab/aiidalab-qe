import pytest


@pytest.fixture()
def workchain_settings_generator():
    """Return a function that generates a workchain settings dictionary."""
    from aiidalab_qe.app.steps import WorkChainSettings

    def _workchain_settings_generator(**kwargs):
        workchain_settings = WorkChainSettings()
        workchain_settings._update_settings(**kwargs)
        return workchain_settings

    return _workchain_settings_generator


@pytest.fixture()
def smearing_settings_generator():
    """Return a function that generates a smearing settings dictionary."""
    from aiidalab_qe.app.steps import SmearingSettings

    def _smearing_settings_generator(**kwargs):
        smearing_settings = SmearingSettings()
        smearing_settings._update_settings(**kwargs)
        print("smearing", smearing_settings.smearing.value)
        return smearing_settings

    return _smearing_settings_generator


@pytest.fixture()
def kpoints_settings_generator():
    """Return a function that generates a kpoints settings dictionary."""
    from aiidalab_qe.app.steps import KpointSettings

    def _kpoints_settings_generator(**kwargs):
        kpoints_settings = KpointSettings()
        kpoints_settings._update_settings(**kwargs)
        return kpoints_settings

    return _kpoints_settings_generator


@pytest.mark.usefixtures("sssp")
def test_create_builder_default(
    data_regression,
    pw_code,
    dos_code,
    projwfc_code,
    structure_data_object,
    workchain_settings_generator,
    smearing_settings_generator,
    kpoints_settings_generator,
):
    """ "Test the creation of the workchain builder.

    metal, non-magnetic
    """
    from aiidalab_qe.app.pseudos import PseudoFamilySelector
    from aiidalab_qe.app.report import _generate_report_html
    from aiidalab_qe.app.steps import SubmitQeAppWorkChainStep

    # XXX: combine the following three lines into one fixture call
    submit_step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    submit_step.input_structure = structure_data_object
    submit_step.pseudo_family_selector = PseudoFamilySelector()

    # XXX: Codes, may also be set in the step constructor
    submit_step.pw_code.value = pw_code.uuid
    submit_step.dos_code.value = dos_code.uuid
    submit_step.projwfc_code.value = projwfc_code.uuid

    # Settings
    submit_step.workchain_settings = workchain_settings_generator(
        relax_type="positions_cell",
        spin_type="none",
        electronic_type="metal",
        bands_run=True,
        pdos_run=True,
        workchain_protocol="moderate",
    )
    submit_step.kpoints_settings = kpoints_settings_generator(kpoints_distance=0.12)
    submit_step.smearing_settings = smearing_settings_generator(
        smearing="methfessel-paxton", degauss=0.015, override_protocol_smearing=True
    )

    builder, extra_parameters = submit_step._create_builder()

    # check and validate the builder
    got = builder_to_readable_dict(builder)

    # XXX: Instead of using data_regression whole builder, the test can be more explicit
    # such as checking all the scf are the same and regressing only the scf.
    data_regression.check(got)

    # test the report can be properly generated from the builder without errors
    builder_parameters = submit_step._extract_report_parameters(
        builder, extra_parameters
    )
    report_html = _generate_report_html(builder_parameters)

    # None in report_html means that the report not properly generated
    assert "None" not in report_html


@pytest.mark.usefixtures("sssp")
def test_create_builder_insulator(
    data_regression,
    pw_code,
    dos_code,
    projwfc_code,
    structure_data_object,
    workchain_settings_generator,
    smearing_settings_generator,
    kpoints_settings_generator,
):
    """ "Test the creation of the workchain builder.

    insulator, non-magnetic, no smearing
    the occupation type is set to fixed, smearing and degauss should not be set"""
    from aiidalab_qe.app.pseudos import PseudoFamilySelector
    from aiidalab_qe.app.report import _generate_report_html
    from aiidalab_qe.app.steps import SubmitQeAppWorkChainStep

    # XXX: combine the following three lines into one fixture call
    submit_step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    submit_step.input_structure = structure_data_object
    submit_step.pseudo_family_selector = PseudoFamilySelector()

    # XXX: Codes, may also be set in the step constructor
    submit_step.pw_code.value = pw_code.uuid
    submit_step.dos_code.value = dos_code.uuid
    submit_step.projwfc_code.value = projwfc_code.uuid

    # Settings
    submit_step.workchain_settings = workchain_settings_generator(
        relax_type="positions_cell",
        spin_type="none",
        electronic_type="insulator",
        bands_run=True,
        pdos_run=True,
        workchain_protocol="moderate",
    )
    submit_step.kpoints_settings = kpoints_settings_generator(kpoints_distance=0.12)
    submit_step.smearing_settings = smearing_settings_generator(
        smearing="methfessel-paxton", degauss=0.015, override_protocol_smearing=True
    )  # override but has not effect

    builder, extra_parameters = submit_step._create_builder()

    # check and validate the builder
    got = builder_to_readable_dict(builder)

    data_regression.check(got)

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
