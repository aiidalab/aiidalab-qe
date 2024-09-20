import traitlets as tl

from aiidalab_widgets_base import WizardAppWidgetStep


class ResultsModel(tl.HasTraits):
    state = tl.UseEnum(
        enum_class=WizardAppWidgetStep.State,
        default_value=WizardAppWidgetStep.State.INIT,
    )
    process = tl.Unicode(allow_none=True)


results_model = ResultsModel()
