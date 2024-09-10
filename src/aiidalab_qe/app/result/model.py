import traitlets

from aiidalab_widgets_base import WizardAppWidgetStep


class ResultsModel(traitlets.HasTraits):
    state = traitlets.UseEnum(WizardAppWidgetStep.State)

    process = traitlets.Unicode(allow_none=True)


results_model = ResultsModel()
