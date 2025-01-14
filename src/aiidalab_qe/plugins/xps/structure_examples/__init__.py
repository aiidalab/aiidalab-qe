import pathlib

file_path = pathlib.Path(__file__).parent
structure_examples = {
    "title": "XPS examples",
    "structures": [
        ("Phenylacetylene molecule", file_path / "Phenylacetylene.xyz"),
        ("ETFA molecule", file_path / "ETFA.xyz"),
    ],
}
