import pathlib

file_path = pathlib.Path(__file__).parent

structure_examples = {
    "title": "XAS examples",
    "structures": [
        ("Lithium carbonate", file_path / "Li2CO3.cif"),
    ],
}
