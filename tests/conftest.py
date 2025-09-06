from __future__ import annotations

import io
import typing as t

import numpy as np
import pytest
from aiida_pseudo.data.pseudo import UpfData
from aiida_pseudo.groups.family import PseudoDojoFamily, SsspFamily

from aiida import engine, orm
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager
from aiida.orm.utils.serialize import serialize
from aiida.plugins import DataFactory, OrbitalFactory
from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiidalab_qe.app.configuration.advanced import (
    AdvancedConfigurationSettingsModel,
    ConvergenceConfigurationSettingsModel,
    GeneralConfigurationSettingsModel,
    MagnetizationConfigurationSettingsModel,
    PseudosConfigurationSettingsModel,
    SmearingConfigurationSettingsModel,
)
from aiidalab_qe.app.configuration.basic import BasicConfigurationSettingsModel
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.wizard import Wizard, WizardModel
from aiidalab_qe.plugins.bands.bands_workchain import BandsWorkChain
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION
from aiidalab_qe.utils import shallow_copy_nested_dict
from aiidalab_qe.workflows import QeAppWorkChain

pytest_plugins = ["aiida.tools.pytest_fixtures"]

ELEMENTS = [
    "H",
    "Li",
    "O",
    "S",
    "Si",
    "Co",
    "Mo",
    "Ni",
]

CUTOFFS = {
    element: {
        "cutoff_wfc": 30.0,
        "cutoff_rho": 240.0,
    }
    for element in ELEMENTS
}


def pytest_addoption(parser):
    parser.addoption(
        "--skip-slow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--skip-slow"):
        return
    skip_slow = pytest.mark.skip(reason="Slow test")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost


@pytest.fixture
def generate_structure_data():
    """generate a `StructureData` object."""

    def _generate_structure_data(name="silicon", pbc=(True, True, True), store=True):
        if name == "silicon":
            structure = orm.StructureData(
                cell=[
                    [3.84737, 0.0, 0.0],
                    [1.923685, 3.331920, 0.0],
                    [1.923685, 1.110640, 3.141364],
                ]
            )
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Si")
            structure.append_atom(position=(1.923685, 1.110640, 0.785341), symbols="Si")

        elif name == "silica":
            structure = orm.StructureData(
                cell=[
                    [4.18, 0.0, 0.0],
                    [0.0, 4.18, 0.0],
                    [0.0, 0.0, 2.66],
                ]
            )
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Si")
            structure.append_atom(position=(2.09, 2.09, 1.33), symbols="Si")
            structure.append_atom(position=(3.37, 0.81, 1.33), symbols="O")
            structure.append_atom(position=(1.28, 1.28, 0.0), symbols="O")
            structure.append_atom(position=(2.9, 2.9, 0.0), symbols="O")
            structure.append_atom(position=(0.81, 3.37, 1.33), symbols="O")

        elif name == "LiCoO2":
            a, b, c, d = (
                1.4060463552647,
                0.81178124180108,
                4.6012019181836,
                1.6235624832021,
            )
            cell = [[a, -b, c], [0.0, d, c], [-a, -b, c]]
            sites = [
                ["Co", "Co", (0, 0, 0)],
                ["O", "O", (0, 0, 3.6020728736387)],
                ["O", "O", (0, 0, 10.201532881212)],
                ["Li", "Li", (0, 0, 6.9018028772754)],
            ]
            structure = orm.StructureData(cell=cell)

            for site in sites:
                structure.append_atom(position=site[2], symbols=site[0], name=site[1])

        elif name == "MoS2":
            cell = [[3.1922, 0, 0], [-1.5961, 2.7646, 0], [0, 0, 13.3783]]
            structure = orm.StructureData(cell=cell)
            structure.append_atom(position=(-0.0, 1.84, 10.03), symbols="Mo")
            structure.append_atom(position=(1.6, 0.92, 8.47), symbols="S")
            structure.append_atom(position=(1.6, 0.92, 11.6), symbols="S")

        elif name == "H2O":
            cell = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
            structure = orm.StructureData(cell=cell)
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols="H")
            structure.append_atom(position=(0.0, 0.0, 1.0), symbols="O")
            structure.append_atom(position=(0.0, 1.0, 0.0), symbols="H")

        elif name == "H2O-larger":
            cell = [[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]]
            structure = orm.StructureData(cell=cell)
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols="H")
            structure.append_atom(position=(0.0, 0.0, 1.0), symbols="O")
            structure.append_atom(position=(0.0, 1.0, 0.0), symbols="H")

        elif name == "CeO":
            cell = [[3.04, 0.0, 1.76], [1.01, 2.87, 1.76], [0.0, 0.0, 3.51]]
            structure = orm.StructureData(cell=cell)
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Ce")
            structure.append_atom(position=(3.51, 2.03, 1.43), symbols="O")

        structure.pbc = pbc

        if store:
            structure.store()

        return structure

    return _generate_structure_data


@pytest.fixture
def generate_xy_data():
    """Return an ``XyData`` instance."""

    def _generate_xy_data(xvals=None, yvals=None, xlabel=None, ylabel=None):
        """Return an ``XyData`` node.
        xvals and yvals are lists, and should have the same length.
        """
        xunits = "n/a"
        yunits = ["n/a"] * len(ylabel)

        xy_node = orm.XyData()
        xy_node.set_x(xvals, xlabel, xunits)
        xy_node.set_y(yvals, ylabel, yunits)
        xy_node.store()
        return xy_node

    return _generate_xy_data


@pytest.fixture
def generate_bands_data():
    """Return a `BandsData` node."""

    def _generate_bands_data():
        """Return a `BandsData` instance with some basic `kpoints` and `bands` arrays."""
        BandsData = DataFactory("core.array.bands")

        kpoints = np.array([[0.0, 0.0, 0.0]])
        bands = np.array([[-5.64024889]])
        bands_data = BandsData()
        bands_data.set_kpoints(kpoints)
        bands_data.set_bands(bands, units="eV")
        bands_data.store()

        return bands_data

    return _generate_bands_data


@pytest.fixture
def generate_projection_data(generate_bands_data):
    """Return an ``ProjectionData`` instance."""

    def _generate_projection_data():
        """Return an ``ProjectionData`` node."""
        ProjectionData = DataFactory("core.array.projection")
        OrbitalCls = OrbitalFactory("core.realhydrogen")

        state_dict = {
            "kind_name": "C",
            "angular_momentum": 0,
            "magnetic_number": 0,
            "radial_nodes": 1,
            "position": [0.0, 0.0, 0.0],
        }
        orbitals = [OrbitalCls(**state_dict)]
        energy_arrays = np.array([1])
        pdos_arrays = np.array([1])

        projection_data = ProjectionData()
        bands_data = generate_bands_data()
        projection_data.set_reference_bandsdata(bands_data)
        projection_data.set_projectiondata(
            orbitals,
            list_of_energy=energy_arrays,
            list_of_pdos=pdos_arrays,
            bands_check=False,
        )
        projection_data.store()
        return projection_data

    return _generate_projection_data


def _generate_upf_data(element, filename=None, params=None):
    """Return `UpfData` node.

    NOTE: `params` is used to render the generated UPF unique per pseudo
    """
    extras = "\n\t\t".join(f'{key}="{value}"' for key, value in (params or {}).items())
    content = f"""
        <UPF version="2.0.1">
            <PP_HEADER
                element="{element}"
                z_valence="4.0"
                {extras}
            />
        </UPF>
    """
    stream = io.BytesIO(content.encode("utf-8"))

    filename = element
    if params:
        filename += f"_{'_'.join(params.values())}"
    filename += ".upf"

    return UpfData(stream, filename=filename)


@pytest.fixture(scope="function")
def generate_upf_data():
    """Return a `UpfData` instance for the given element a file for which should exist in `tests/fixtures/pseudos`."""
    return _generate_upf_data


@pytest.fixture(scope="session")
def generate_upf_data_for_session():
    """Return a `UpfData` instance for the given element a file for which should exist in `tests/fixtures/pseudos`.

    NOTE: Duplicated fixture for use in the session-scoped pseudos fixtures
    """
    return _generate_upf_data


@pytest.fixture(scope="session", autouse=True)
def sssp(generate_upf_data_for_session):
    """Create SSSP pseudopotentials from scratch."""

    for functional in ("PBE", "PBEsol"):
        for accuracy in ("efficiency", "precision"):
            label = f"SSSP/{SSSP_VERSION}/{functional}/{accuracy}"
            family = SsspFamily(label)
            family.store()
            nodes = []
            for element in ELEMENTS:
                params = {
                    "functional": functional,
                    "accuracy": accuracy,
                }
                filename = f"{element}_{'_'.join(params.values())}.upf"
                upf = generate_upf_data_for_session(element, filename, params)
                upf.store()
                nodes.append(upf)
            family.add_nodes(nodes)
            family.set_cutoffs(CUTOFFS, accuracy, unit="Ry")


@pytest.fixture(scope="session", autouse=True)
def pseudodojo(generate_upf_data_for_session):
    """Create pseudodojo pseudopotentials from scratch."""

    ROOT = f"PseudoDojo/{PSEUDODOJO_VERSION}"
    for stringency in ("standard", "stringent"):
        for functional in ("PBE", "PBEsol"):
            for relativistic in ("SR", "FR"):
                label = f"{ROOT}/{functional}/{relativistic}/{stringency}/upf"
                family = PseudoDojoFamily(label)
                family.store()
                nodes = []
                for element in ELEMENTS:
                    params = {
                        "functional": functional,
                        "stringency": stringency,
                        "relativistic": relativistic,
                    }
                    filename = f"{element}_{'_'.join(params.values())}.upf"
                    upf = generate_upf_data_for_session(element, filename, params)
                    upf.store()
                    nodes.append(upf)
                family.add_nodes(nodes)
                family.set_cutoffs(CUTOFFS, stringency, unit="Ry")


@pytest.fixture
def generate_code(aiida_code_installed):
    def _generate_code(label):
        return aiida_code_installed(
            label=label,
            default_calc_job_plugin=f"quantumespresso.{label}",
        )

    return _generate_code


@pytest.fixture
def pw_code(generate_code):
    return generate_code("pw")


@pytest.fixture
def dos_code(generate_code):
    return generate_code("dos")


@pytest.fixture
def projwfc_code(generate_code):
    return generate_code("projwfc")


@pytest.fixture
def app(pw_code, dos_code, projwfc_code):
    # Assign test codes as defaults
    DEFAULTS = t.cast(dict, DEFAULT_PARAMETERS)
    DEFAULTS["codes"]["pw"]["code"] = pw_code.full_label
    DEFAULTS["codes"]["dos"]["code"] = dos_code.full_label
    DEFAULTS["codes"]["projwfc"]["code"] = projwfc_code.full_label

    model = WizardModel()
    app = Wizard(model=model, auto_setup=False)

    # Since we use `auto_setup=False`, which will skip the pseudo library
    # installation, we need to mock set the installation status to `True` to
    # avoid the blocker message pop up in the submission step.
    app.structure_model.installing_sssp = False
    app.structure_model.sssp_installed = True
    app.submit_model.installing_qe = False
    app.submit_model.qe_installed = True

    yield app


@pytest.fixture()
def submit_app_generator(
    app: Wizard,
    generate_structure_data,
):
    """Return a function that generates a submit step widget."""

    def _submit_app_generator(
        relax_type="positions_cell",
        spin_type="none",
        electronic_type="metal",
        properties=None,
        workchain_protocol="balanced",
        kpoints_distance=0.12,
        smearing="methfessel-paxton",
        degauss=0.015,
        tot_charge=0.0,
        vdw_corr="none",
        initial_magnetic_moments=0.0,
        electron_maxstep=80,
    ):
        app.structure_model.structure_uuid = generate_structure_data().uuid
        app.structure_model.confirm()

        parameters = {
            "workchain": {
                "relax_type": relax_type,
                "spin_type": spin_type,
                "electronic_type": electronic_type,
                "properties": properties or [],
                "protocol": workchain_protocol,
            }
        }
        app.configure_model.set_model_state(parameters)

        advanced_model = t.cast(
            AdvancedConfigurationSettingsModel,
            app.configure_model.get_model("advanced"),
        )

        general_model = t.cast(
            GeneralConfigurationSettingsModel,
            advanced_model.get_model("general"),
        )
        general_model.total_charge = tot_charge
        general_model.van_der_waals = vdw_corr

        convergence_model = t.cast(
            ConvergenceConfigurationSettingsModel,
            advanced_model.get_model("convergence"),
        )
        convergence_model.electron_maxstep = electron_maxstep
        convergence_model.kpoints_distance = kpoints_distance

        smearing_model = t.cast(
            SmearingConfigurationSettingsModel,
            advanced_model.get_model("smearing"),
        )
        smearing_model.type = smearing
        smearing_model.degauss = degauss

        if isinstance(initial_magnetic_moments, (int, float)):
            initial_magnetic_moments = [initial_magnetic_moments]

        magnetization_model = t.cast(
            MagnetizationConfigurationSettingsModel,
            advanced_model.get_model("magnetization"),
        )
        magnetization_model.moments = dict(
            zip(
                app.configure_model.input_structure.get_kind_names(),
                initial_magnetic_moments,
            )
        )

        app.configure_model.confirm()

        global_resources_model = app.submit_model.get_model("global")
        global_resources_model.get_model("quantumespresso__pw").num_cpus = 2

        return app

    return _submit_app_generator


@pytest.fixture
def app_to_submit(app: Wizard, generate_structure_data):
    # Step 1: select structure from example
    app.structure_model.structure_uuid = generate_structure_data().uuid
    app.structure_model.confirm()
    # Step 2: configure calculation
    # TODO do we need to include bands and pdos here?
    app.configure_model.get_model("bands").include = True
    app.configure_model.get_model("pdos").include = True
    app.configure_model.confirm()
    yield app


@pytest.fixture
def generate_workchain():
    """Generate an instance of a `WorkChain`."""

    def _generate_workchain(process_class, inputs):
        """Generate an instance of a `WorkChain` with the given entry point and inputs.

        :param entry_point: entry point name of the work chain subclass.
        :param inputs: inputs to be passed to process construction.
        :return: a `WorkChain` instance.
        """
        runner = get_manager().get_runner()
        process = instantiate_process(runner, process_class, **inputs)

        return process

    return _generate_workchain


@pytest.fixture
def generate_pdos_workchain(
    fixture_localhost,
    pw_code,
    dos_code,
    projwfc_code,
    generate_xy_data,
    generate_projection_data,
    generate_workchain,
):
    """Generate an instance of a `XpsWorkChain`."""

    def _generate_pdos_workchain(structure, spin_type="none"):
        pseudo_family = f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"

        inputs = {
            "pw_code": pw_code,
            "dos_code": dos_code,
            "projwfc_code": projwfc_code,
            "structure": structure,
            "overrides": {
                "scf": {
                    "pseudo_family": pseudo_family,
                },
                "nscf": {
                    "pseudo_family": pseudo_family,
                },
            },
        }
        builder = PdosWorkChain.get_builder_from_protocol(**inputs)
        inputs = builder._inputs()
        wkchain = generate_workchain(PdosWorkChain, inputs)
        wkchain.setup()
        # wkchain.run_pdos()
        remote = orm.RemoteData(remote_path="/tmp/aiida_run")
        remote.computer = fixture_localhost
        remote.store()
        retrieved = orm.FolderData(tree="/tmp/aiida_run")
        retrieved.store()
        output_parameters = orm.Dict(dict={"fermi_energy": 2.0})
        output_parameters.store()
        proj = generate_projection_data()
        proj.store()
        if spin_type == "none":
            xy = generate_xy_data(
                np.array([1, 2, 3]), [np.array([1, 2, 3])], "X", ["dos"]
            )
            xy.store()
            wkchain.out(
                "dos",
                {
                    "output_dos": xy,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
            wkchain.out(
                "projwfc",
                {
                    "Dos": xy,
                    "projections": proj,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
        else:
            xy = generate_xy_data(
                np.array([1, 2, 3]),
                [np.array([1, 2, 3]), np.array([1, 2, 3])],
                "X",
                ["dos_spin_up", "dos_spin_down"],
            )
            xy.store()
            wkchain.out(
                "dos",
                {
                    "output_dos": xy,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
            wkchain.out(
                "projwfc",
                {
                    "Dos": xy,
                    "projections_up": proj,
                    "projections_down": proj,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
        wkchain.out(
            "nscf",
            {
                "output_parameters": output_parameters,
                "remote_folder": remote,
                "retrieved": retrieved,
            },
        )
        wkchain.update_outputs()
        pdos_node = wkchain.node
        pdos_node.set_exit_status(0)
        pdos_node.set_process_state(engine.ProcessState.FINISHED)
        # set
        return wkchain

    return _generate_pdos_workchain


@pytest.fixture
def generate_bands_workchain(
    pw_code,
    projwfc_code,
    generate_bands_data,
    generate_workchain,
):
    """Generate an instance of a the WorkChain."""

    def _generate_bands_workchain(structure):
        pseudo_family = f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"

        inputs = {
            "pw_code": pw_code,
            "projwfc_code": projwfc_code,
            "structure": structure,
            "simulation_mode": "normal",
            "overrides": {
                "scf": {
                    "pseudo_family": pseudo_family,
                },
                "bands": {
                    "pseudo_family": pseudo_family,
                },
                "relax": {
                    "base": {
                        "pseudo_family": pseudo_family,
                    },
                    "base_final_scf": {
                        "pseudo_family": pseudo_family,
                    },
                },
            },
        }
        builder = BandsWorkChain.get_builder_from_protocol(**inputs)
        inputs = builder._inputs()
        wkchain = generate_workchain(BandsWorkChain, inputs)
        wkchain.setup()
        # run bands and return the process
        fermi_dict = orm.Dict(dict={"fermi_energy": 2.0})
        fermi_dict.store()
        output_parameters = {
            "bands": {
                "scf_parameters": fermi_dict,
                "band_parameters": fermi_dict,
            }
        }

        wkchain.out(
            "bands.scf_parameters", output_parameters["bands"]["scf_parameters"]
        )
        wkchain.out(
            "bands.band_parameters", output_parameters["bands"]["band_parameters"]
        )

        #
        band_structure = generate_bands_data()
        band_structure.store()
        wkchain.out("bands.band_structure", band_structure)
        wkchain.update_outputs()
        #
        bands_node = wkchain.node
        bands_node.set_exit_status(0)
        bands_node.set_process_state(engine.ProcessState.FINISHED)
        return wkchain

    return _generate_bands_workchain


@pytest.fixture
def generate_qeapp_workchain(
    app: Wizard,
    projwfc_code,
    generate_structure_data,
    generate_workchain,
    generate_pdos_workchain,
    generate_bands_workchain,
):
    """Generate an instance of the WorkChain."""

    def _generate_qeapp_workchain(
        structure: orm.StructureData | None = None,
        relax_type="positions_cell",
        run_bands=True,
        run_pdos=True,
        spin_type="none",
        electronic_type="metal",
        magnetization_type="starting_magnetization",  # Options: "starting_magnetization", "tot_magnetization"
        initial_magnetic_moments=0.0,
        tot_magnetization=0.0,
        functional="PBEsol",
    ):
        # Step 1: select structure from example
        if structure is None:
            structure = generate_structure_data()
        else:
            structure.store()

        app.structure_model.structure_uuid = structure.uuid
        app.structure_model.confirm()

        # step 2 configure
        basic_model = t.cast(
            BasicConfigurationSettingsModel,
            app.configure_model.get_model("workchain"),
        )

        app.configure_model.relax_type = relax_type

        # In order to prepare complete inputs, I set all the properties to true
        # this can be overridden later
        app.configure_model.get_model("bands").include = run_bands
        app.configure_model.get_model("pdos").include = run_pdos

        basic_model.protocol = "fast"
        basic_model.spin_type = spin_type
        basic_model.electronic_type = electronic_type

        advanced_model = t.cast(
            AdvancedConfigurationSettingsModel,
            app.configure_model.get_model("advanced"),
        )

        if spin_type == "collinear":
            magnetization_model = t.cast(
                MagnetizationConfigurationSettingsModel,
                advanced_model.get_model("magnetization"),
            )
            if electronic_type == "insulator":
                magnetization_model.total = tot_magnetization
            elif magnetization_type == "starting_magnetization":
                if isinstance(initial_magnetic_moments, (int, float)):
                    initial_magnetic_moments = [initial_magnetic_moments]
                magnetization_model.moments = dict(
                    zip(
                        structure.get_kind_names(),
                        initial_magnetic_moments,
                    )
                )
            else:
                magnetization_model.total = tot_magnetization

        pseudos = t.cast(
            PseudosConfigurationSettingsModel,
            advanced_model.get_model("pseudos"),
        )
        pseudos.functional = functional

        app.configure_model.confirm()

        # step 3 setup code and resources
        global_resources_model = app.submit_model.get_model("global")
        global_resources_model.get_model("quantumespresso__pw").num_cpus = 4
        parameters = app.submit_model.get_model_state()
        builder = app.submit_model._create_builder(parameters)

        inputs = builder._inputs()
        if "relax" in inputs:
            inputs["relax"]["base_final_scf"] = shallow_copy_nested_dict(
                inputs["relax"]["base"]
            )

        if run_bands:
            # Setting up inputs for bands_projwfc
            inputs["bands"]["bands_projwfc"]["scf"]["pw"] = shallow_copy_nested_dict(
                inputs["bands"]["bands"]["scf"]["pw"]
            )
            inputs["bands"]["bands_projwfc"]["bands"]["pw"] = shallow_copy_nested_dict(
                inputs["bands"]["bands"]["bands"]["pw"]
            )
            inputs["bands"]["bands_projwfc"]["bands"]["pw"]["code"] = inputs["bands"][
                "bands"
            ]["bands"]["pw"]["code"]
            inputs["bands"]["bands_projwfc"]["scf"]["pw"]["code"] = inputs["bands"][
                "bands"
            ]["scf"]["pw"]["code"]

            inputs["bands"]["bands_projwfc"]["projwfc"]["projwfc"]["code"] = (
                projwfc_code
            )
            inputs["bands"]["bands_projwfc"]["projwfc"]["projwfc"]["parameters"] = (
                orm.Dict({"PROJWFC": {"DeltaE": 0.01}}).store()
            )
            inputs["properties"].append("bands")

        if run_pdos:
            inputs["properties"].append("pdos")

        workchain = generate_workchain(QeAppWorkChain, inputs)
        workchain.setup()

        # mock output
        if relax_type != "none":
            workchain.out("structure", app.structure_model.input_structure)

        if run_pdos:
            pdos = generate_pdos_workchain(structure, spin_type)
            workchain.out_many(
                workchain.exposed_outputs(
                    pdos.node,
                    PdosWorkChain,
                    namespace="pdos",
                )
            )

        if run_bands:
            bands = generate_bands_workchain(structure)
            workchain.out_many(
                workchain.exposed_outputs(
                    bands.node,
                    BandsWorkChain,
                    namespace="bands",
                )
            )
        workchain.update_outputs()

        # set ui_parameters
        qeapp_node = workchain.node
        qeapp_node.base.extras.set("ui_parameters", serialize(parameters))

        return workchain

    return _generate_qeapp_workchain
