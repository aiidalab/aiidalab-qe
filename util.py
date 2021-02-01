import json

from aiida.orm import Dict


DEFAULT_PARAMETERS_FILE = "parameters.json"


def load_default_parameters():
    with open(DEFAULT_PARAMETERS_FILE) as file:
        return Dict(dict=json.load(file))
