import traitlets as tl


class ResultsModel(tl.HasTraits):
    process = tl.Unicode(allow_none=True)

    process_info = tl.Unicode("")

    def reset(self):
        self.process = None
        self.process_info = ""
