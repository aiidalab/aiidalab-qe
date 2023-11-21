import tempfile

from aiida.orm import CalcJobNode
from aiida.transports import Transport


def relax_pw(node: CalcJobNode, transport: Transport) -> str:
    """Retrieve and save the total energy from aiida.out file and save it into the node extras.

    :param node: The node representing the calculation job.
    :param transport: The transport that can be used to retrieve files from remote working directory.
    :returns: A string if the job should be killed, `None` otherwise.
    """
    import numpy as np
    from ase.io import read

    try:
        # read structure
        with tempfile.NamedTemporaryFile("w+") as handle:
            transport.getfile(node.get_option("output_filename"), handle.name)
            handle.seek(0)
            images = read(handle, index=":", format="espresso-out")
            if len(images) > 0:
                atoms = images[0].todict()
                # type `<class 'numpy.bool_'>` is not supported as it is not json-serializable
                atoms.pop("pbc")
                trajectory = np.zeros((len(images), len(images[0]), 3))
                for i in range(len(images)):
                    trajectory[i] = images[i].positions
                node.base.extras.set("monitor_atoms", atoms)
                node.base.extras.set("monitor_trajectory", trajectory)
        # read scf energy
        with tempfile.NamedTemporaryFile("w+") as handle:
            transport.getfile(node.get_option("output_filename"), handle.name)
            handle.seek(0)
            output = handle.read()
        energies = []
        lines = output.splitlines()
        for line in lines:
            if "total energy              =" in line:
                energies.append(float(line.split("=")[1].split()[0]))
        node.base.extras.set("monitor_energies", energies)
    except FileExistsError:
        return
