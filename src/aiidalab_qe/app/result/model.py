import traitlets as tl


class ResultsModel(tl.HasTraits):
    process = tl.Unicode(allow_none=True)


results_model = ResultsModel()
