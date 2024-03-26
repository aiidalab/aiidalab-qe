from aiida_quantumespresso_hp.workflows.hp.main import HpWorkChain


def check_codes(pw_code, hp_code):
    """Check that the codes are installed on the same computer."""
    if (
        not any(
            [
                pw_code is None,
                hp_code is None,
            ]
        )
        and len(
            set(
                (
                    pw_code.computer.pk,
                    hp_code.computer.pk,
                )
            )
        )
        != 1
    ):
        raise ValueError(
            "All selected codes must be installed on the same computer. This is because the "
            "HP calculations rely on large files that are not retrieved by AiiDA."
        )


def get_builder(codes, structure, parameters, **kwargs):
    pw_code = codes.get("pw")
    hp_code = codes.get("hp")
    check_codes(pw_code, hp_code)
    protocol = parameters["workchain"]["protocol"]

    overrides = {
        "parallelize_atoms": True,
        "parallelize_qpoints": True,
        "hp": {
            "hubbard_structure": structure,
            "metadata": {
                "options": {
                    "resources": {
                        "num_machines": 1,
                        "num_mpiprocs_per_machine": 2,
                    },
                }
            },
        },
        "qpoints_distance": 1000,  # to get few q points
    }
    if hp_code is not None:
        builder = HpWorkChain.get_builder_from_protocol(
            code=hp_code,
            protocol=protocol,
            overrides=overrides,
            **kwargs,
        )
    else:
        raise ValueError("The hp_code and bader_code are required.")
    return builder


workchain_and_builder = {
    "workchain": HpWorkChain,
    "exclude": ("clean_workdir", "structure", "relax"),
    "get_builder": get_builder,
    "input_from_ctx": {"hp.parent_folder": "scf_folder"},
}
