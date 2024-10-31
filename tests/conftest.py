from __future__ import annotations

import io
import pathlib
import tempfile

import pytest

from aiida import orm
from aiidalab_qe.app.configuration.model import ConfigurationModel
from aiidalab_qe.app.main import App
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION

pytest_plugins = ["aiida.manage.tests.pytest_fixtures"]

ELEMENTS = [
    "H",
    "Li",
    "O",
    "S",
    "Si",
    "Co",
    "Mo",
]


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
def fixture_code(fixture_localhost):
    """Return an ``InstalledCode`` instance configured to run calculations of given entry point on localhost."""

    def _fixture_code(entry_point_name):
        from aiida.orm import InstalledCode, load_code

        label = f"test.{entry_point_name}"

        try:
            return load_code(label=label)
        except Exception:
            return InstalledCode(
                label=label,
                computer=fixture_localhost,
                filepath_executable="/bin/true",
                default_calc_job_plugin=entry_point_name,
            )

    return _fixture_code


@pytest.fixture
def generate_structure_data():
    """generate a `StructureData` object."""

    def _generate_structure_data(name="silicon", pbc=(True, True, True)):
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

        structure.pbc = pbc

        return structure

    return _generate_structure_data


@pytest.fixture
def generate_xy_data():
    """Return an ``XyData`` instance."""

    def _generate_xy_data(xvals=None, yvals=None, xlabel=None, ylabel=None):
        """Return an ``XyData`` node.
        xvals and yvals are lists, and should have the same length.
        """
        from aiida.orm import XyData

        xunits = "n/a"
        yunits = ["n/a"] * len(ylabel)

        xy_node = XyData()
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
        import numpy as np

        from aiida.plugins import DataFactory

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
        import numpy as np

        from aiida.plugins import DataFactory, OrbitalFactory

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


@pytest.fixture(scope="function")
def sssp(generate_upf_data):
    """Create SSSP pseudopotentials from scratch."""
    from aiida_pseudo.groups.family import SsspFamily

    cutoffs = {}
    stringency = "standard"

    with tempfile.TemporaryDirectory() as d:
        dirpath = pathlib.Path(d)

        for element in ELEMENTS:
            upf = generate_upf_data(element)
            filename = dirpath / f"{element}.upf"

            with open(filename, "w+b") as handle:
                with upf.open(mode="rb") as source:
                    handle.write(source.read())
                    handle.flush()

            cutoffs[element] = {
                "cutoff_wfc": 30.0,
                "cutoff_rho": 240.0,
            }

        for functional in ("PBE", "PBEsol"):
            for accuracy in ("efficiency", "precision"):
                label = f"SSSP/{SSSP_VERSION}/{functional}/{accuracy}"
                family = SsspFamily.create_from_folder(dirpath, label)
                family.set_cutoffs(cutoffs, stringency, unit="Ry")


@pytest.fixture(scope="function")
def pseudodojo(generate_upf_data):
    """Create pseudodojo pseudopotentials from scratch."""
    from aiida_pseudo.data.pseudo import UpfData
    from aiida_pseudo.groups.family import PseudoDojoFamily

    cutoffs = {}

    with tempfile.TemporaryDirectory() as d:
        dirpath = pathlib.Path(d)

        for element in ELEMENTS:
            upf = generate_upf_data(element)
            filename = dirpath / f"{element}.upf"

            upf = generate_upf_data(element)
            filename = dirpath / f"{element}.upf"

            with open(filename, "w+b") as handle:
                with upf.open(mode="rb") as source:
                    handle.write(source.read())
                    handle.flush()

            cutoffs[element] = {
                "cutoff_wfc": 36.0,
                "cutoff_rho": 144.0,
            }

        ROOT = f"PseudoDojo/{PSEUDODOJO_VERSION}"
        for stringency in ("standard", "stringent"):
            for functional in ("PBE", "PBEsol"):
                for spin_orbit in ("SR", "FR"):
                    label = f"{ROOT}/{functional}/{spin_orbit}/{stringency}/upf"
                    family = PseudoDojoFamily.create_from_folder(
                        dirpath,
                        label,
                        pseudo_type=UpfData,
                    )
                    family.set_cutoffs(cutoffs, stringency, unit="Ry")


@pytest.fixture(scope="function")
def generate_upf_data():
    """Return a `UpfData` instance for the given element a file for which should exist in `tests/fixtures/pseudos`."""

    def _generate_upf_data(element, filename=None):
        """Return `UpfData` node."""
        from aiida_pseudo.data.pseudo import UpfData

        content = f'<UPF version="2.0.1"><PP_HEADER\nelement="{element}"\nz_valence="4.0"\n/></UPF>\n'
        stream = io.BytesIO(content.encode("utf-8"))

        if filename is None:
            filename = f"{element}.upf"

        return UpfData(stream, filename=filename)

    return _generate_upf_data


@pytest.fixture
def pw_code(aiida_local_code_factory):
    """Return a `Code` configured for the pw.x executable."""

    return aiida_local_code_factory(
        label="pw", executable="bash", entry_point="quantumespresso.pw"
    )


@pytest.fixture
def dos_code(aiida_local_code_factory):
    """Return a `Code` configured for the dos.x executable."""
    return aiida_local_code_factory(
        label="dos", executable="bash", entry_point="quantumespresso.dos"
    )


@pytest.fixture
def projwfc_code(aiida_local_code_factory):
    """Return a `Code` configured for the projwfc.x executable."""
    return aiida_local_code_factory(
        label="projwfc",
        executable="bash",
        entry_point="quantumespresso.projwfc",
    )


@pytest.fixture
def projwfc_bands_code(aiida_local_code_factory):
    """Return a `Code` configured for the projwfc.x executable."""
    return aiida_local_code_factory(
        label="projwfc_bands",
        executable="bash",
        entry_point="quantumespresso.projwfc",
    )


@pytest.fixture()
def workchain_settings_generator():
    """Return a function that generates a workchain settings dictionary."""
    from aiidalab_qe.app.configuration.basic.workflow import BasicSettings

    def _workchain_settings_generator(**kwargs):
        model = ConfigurationModel()
        workchain_settings = BasicSettings(config_model=model)
        workchain_settings._update_settings(**kwargs)
        return workchain_settings

    return _workchain_settings_generator


@pytest.fixture()
def smearing_settings_generator():
    """Return a function that generates a smearing settings dictionary."""
    from aiidalab_qe.app.configuration.advanced.smearing import SmearingSettings

    def _smearing_settings_generator(**kwargs):
        model = ConfigurationModel()
        smearing_settings = SmearingSettings(model=model)
        smearing_settings.update_settings(**kwargs)
        return smearing_settings

    return _smearing_settings_generator


@pytest.fixture
def app(pw_code, dos_code, projwfc_code, projwfc_bands_code):
    from aiidalab_qe.app.main import App

    app = App(qe_auto_setup=False)
    app.structure_step.render()
    app.configure_step.render()
    app.submit_step.render()
    app.results_step.render()

    # Since we use `qe_auto_setup=False`, which will skip the pseudo library
    # installation, we need to mock set the installation status to `True` to
    # avoid the blocker message pop up in the submission step.
    app.submit_step.sssp_installation.installed = True
    app.submit_step.qe_setup.installed = True

    # set up codes
    app.submit_model.get_code("pdos", "dos").activate()
    app.submit_model.get_code("pdos", "projwfc").activate()
    app.submit_model.get_code("bands", "projwfc_bands").activate()

    app.submit_model.code_widgets["pw"].code_selection.refresh()
    app.submit_model.code_widgets["dos"].code_selection.refresh()
    app.submit_model.code_widgets["projwfc"].code_selection.refresh()
    app.submit_model.code_widgets["projwfc_bands"].code_selection.refresh()

    app.submit_model.code_widgets["pw"].value = pw_code.uuid
    app.submit_model.code_widgets["dos"].value = dos_code.uuid
    app.submit_model.code_widgets["projwfc"].value = projwfc_code.uuid
    app.submit_model.code_widgets["projwfc_bands"].value = projwfc_bands_code.uuid

    # TODO overrides app defaults - check!
    app.submit_model.set_selected_codes(
        {
            "pw": {"code": pw_code.label},
            "dos": {"code": dos_code.label},
            "projwfc": {"code": projwfc_code.label},
            "projwfc_bands": {"code": projwfc_bands_code.label},
        }
    )

    yield app


@pytest.fixture()
def submit_app_generator(
    app: App,
    generate_structure_data,
):
    """Return a function that generates a submit step widget."""

    def _submit_app_generator(
        relax_type="positions_cell",
        spin_type="none",
        electronic_type="metal",
        properties=None,
        workchain_protocol="moderate",
        kpoints_distance=0.12,
        smearing="methfessel-paxton",
        degauss=0.015,
        tot_charge=0.0,
        vdw_corr="none",
        initial_magnetic_moments=0.0,
        electron_maxstep=80,
    ):
        # Settings
        app.configure_model.input_structure = generate_structure_data()
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

        # Advanced settings
        advanced_model = app.configure_model.get_model("advanced")

        advanced_model.override = True
        advanced_model.total_charge = tot_charge
        advanced_model.van_der_waals = vdw_corr
        advanced_model.kpoints_distance = kpoints_distance
        advanced_model.electron_maxstep = electron_maxstep
        if isinstance(initial_magnetic_moments, (int, float)):
            initial_magnetic_moments = [initial_magnetic_moments]
        advanced_model.get_model("magnetization").moments = dict(
            zip(
                app.configure_model.input_structure.get_kind_names(),
                initial_magnetic_moments,
            )
        )
        # mimic the behavior of the smearing widget set up
        smearing_model = advanced_model.get_model("smearing")
        smearing_model.type = smearing
        smearing_model.degauss = degauss
        app.configure_step.confirm()

        app.submit_model.input_structure = generate_structure_data()
        app.submit_model.code_widgets["pw"].num_cpus.value = 2

        return app

    return _submit_app_generator


@pytest.fixture
def app_to_submit(app: App):
    # Step 1: select structure from example
    structure = app.structure_step.manager.children[0].children[3]  # type: ignore
    structure.children[0].value = structure.children[0].options[1][1]
    app.structure_step.confirm()
    # Step 2: configure calculation
    # TODO do we need to include bands and pdos here?
    app.configure_model.get_model("bands").include = True
    app.configure_model.get_model("pdos").include = True
    app.configure_step.confirm()
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
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager

        runner = get_manager().get_runner()
        process = instantiate_process(runner, process_class, **inputs)

        return process

    return _generate_workchain


@pytest.fixture
def generate_pdos_workchain(
    fixture_localhost,
    fixture_code,
    generate_xy_data,
    generate_projection_data,
    generate_workchain,
):
    """Generate an instance of a `XpsWorkChain`."""

    def _generate_pdos_workchain(structure, spin_type="none"):
        import numpy as np

        from aiida import engine
        from aiida.orm import Dict, FolderData, RemoteData
        from aiida_quantumespresso.workflows.pdos import PdosWorkChain

        pseudo_family = f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"

        inputs = {
            "pw_code": fixture_code("quantumespresso.pw"),
            "dos_code": fixture_code("quantumespresso.dos"),
            "projwfc_code": fixture_code("quantumespresso.projwfc"),
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
        remote = RemoteData(remote_path="/tmp/aiida_run")
        remote.computer = fixture_localhost
        remote.store()
        retrieved = FolderData(tree="/tmp/aiida_run")
        retrieved.store()
        output_parameters = Dict(dict={"fermi_energy": 2.0})
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
    fixture_code,
    generate_bands_data,
    generate_workchain,
):
    """Generate an instance of a the WorkChain."""

    def _generate_bands_workchain(structure):
        from aiida import engine
        from aiida.orm import Dict
        from aiidalab_qe.plugins.bands.bands_workchain import BandsWorkChain

        pseudo_family = f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"

        inputs = {
            "pw_code": fixture_code("quantumespresso.pw"),
            "projwfc_code": fixture_code("quantumespresso.projwfc"),
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
        fermi_dict = Dict(dict={"fermi_energy": 2.0})
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
    app: App,
    generate_workchain,
    generate_pdos_workchain,
    generate_bands_workchain,
    fixture_code,
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
    ):
        from copy import deepcopy

        from aiida.orm import Dict
        from aiida.orm.utils.serialize import serialize
        from aiidalab_qe.workflows import QeAppWorkChain

        # Step 1: select structure from example
        if structure is None:
            from_example = app.structure_step.manager.children[0].children[3]  # type: ignore
            # TODO: (unkpcz) using options to set value in test is cranky, instead, use fixture which will make the test more static and robust.
            from_example.children[0].value = from_example.children[0].options[1][1]
        else:
            structure.store()
            aiida_database_wrapper = app.structure_step.manager.children[0].children[2]  # type: ignore
            aiida_database_wrapper.render()
            aiida_database = aiida_database_wrapper.children[0]  # type: ignore
            aiida_database.search()
            aiida_database.results.value = structure

        app.structure_step.confirm()

        structure = app.structure_model.structure  # type: ignore

        # step 2 configure
        workchain_model = app.configure_model.get_model("workchain")
        advanced_model = app.configure_model.get_model("advanced")

        app.configure_model.relax_type = relax_type

        # In order to prepare complete inputs, I set all the properties to true
        # this can be overridden later
        app.configure_model.get_model("bands").include = run_bands
        app.configure_model.get_model("pdos").include = run_pdos

        workchain_model.protocol = "fast"
        workchain_model.spin_type = spin_type
        workchain_model.electronic_type = electronic_type

        if spin_type == "collinear":
            advanced_model.override = True
            magnetization_model = advanced_model.get_model("magnetization")
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

        app.configure_step.confirm()

        # step 3 setup code and resources
        app.submit_model.code_widgets["pw"].num_cpus.value = 4
        parameters = app.submit_model._get_submission_parameters()
        builder = app.submit_model._create_builder(parameters)

        inputs = builder._inputs()
        inputs["relax"]["base_final_scf"] = deepcopy(inputs["relax"]["base"])

        # Setting up inputs for bands_projwfc
        inputs["bands"]["bands_projwfc"]["scf"]["pw"] = deepcopy(
            inputs["bands"]["bands"]["scf"]["pw"]
        )
        inputs["bands"]["bands_projwfc"]["bands"]["pw"] = deepcopy(
            inputs["bands"]["bands"]["bands"]["pw"]
        )
        inputs["bands"]["bands_projwfc"]["bands"]["pw"]["code"] = inputs["bands"][
            "bands"
        ]["bands"]["pw"]["code"]
        inputs["bands"]["bands_projwfc"]["scf"]["pw"]["code"] = inputs["bands"][
            "bands"
        ]["scf"]["pw"]["code"]

        inputs["bands"]["bands_projwfc"]["projwfc"]["projwfc"]["code"] = fixture_code(
            "quantumespresso.projwfc"
        )
        inputs["bands"]["bands_projwfc"]["projwfc"]["projwfc"]["parameters"] = Dict(
            {"PROJWFC": {"DeltaE": 0.01}}
        ).store()

        if run_bands:
            inputs["properties"].append("bands")
        if run_pdos:
            inputs["properties"].append("pdos")

        workchain = generate_workchain(QeAppWorkChain, inputs)
        workchain.setup()

        # mock output
        if relax_type != "none":
            workchain.out("structure", app.structure_model.structure)
        if run_pdos:
            from aiida_quantumespresso.workflows.pdos import PdosWorkChain

            pdos = generate_pdos_workchain(structure, spin_type)
            workchain.out_many(
                workchain.exposed_outputs(
                    pdos.node,
                    PdosWorkChain,
                    namespace="pdos",
                )
            )
        if run_bands:
            from aiidalab_qe.plugins.bands.bands_workchain import BandsWorkChain

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
