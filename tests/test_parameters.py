from aiidalab_qe.parameters import _set_default_code_labels
from aiidalab_qe.setup.codes import QE_VERSION


def test_set_default_code_labels_injects_qe_version():
    default_codes = {
        "pw": {"code": "pw-{QE_VERSION}@localhost"},
        "projwfc": {"code": "projwfc-{QE_VERSION}@localhost"},
    }

    assert _set_default_code_labels(default_codes) == {
        "pw": {"code": f"pw-{QE_VERSION}@localhost"},
        "projwfc": {"code": f"projwfc-{QE_VERSION}@localhost"},
    }
