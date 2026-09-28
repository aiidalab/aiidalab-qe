import ipywidgets as ipw
import pytest
import traitlets

from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.common.widgets import ProgressBar


def test_progress_bar_defaults_and_assets():
    progress = ProgressBar()

    assert progress.value == 0
    assert progress.bar_style == ""
    assert progress.animating is False
    assert "export const render" in progress._esm
    assert ".qe-progress-bar" in progress._css


def test_progress_bar_preserves_description_layout():
    progress = ProgressBar(description_layout=ipw.Layout(min_width="300px"))

    assert progress.description_layout == {"min_width": "300px"}


@pytest.mark.parametrize("value", [-0.1, 1.1])
def test_progress_bar_rejects_out_of_range_values(value):
    progress = ProgressBar()

    with pytest.raises(traitlets.TraitError, match=r"between 0 and 1\.0"):
        progress.value = value


def test_progress_bar_traits_switch_between_determinate_and_indeterminate():
    progress = ProgressBar()
    states = []
    progress.observe(lambda change: states.append(change["new"]), "animating")

    progress.description = "Installing"
    progress.bar_style = "info"
    progress.animating = True
    progress.value = 0.4
    progress.animating = False

    assert progress.description == "Installing"
    assert progress.bar_style == "info"
    assert progress.value == 0.4
    assert states == [True, False]


def test_qe_setup_widget_stops_animation_when_setup_completes():
    widget = QESetupWidget(auto_start=False)

    widget.set_trait("busy", True)
    assert widget._progress_bar.animating is True

    widget.set_trait("installed", True)
    widget.set_trait("busy", False)
    assert widget._progress_bar.animating is False
    assert widget._progress_bar.value == 1.0
