import ipywidgets as ipw
import traitlets as tl


class AdvancedSubModel(tl.HasTraits):
    dependencies = []

    _defaults = {}

    def update(self, which):
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

    def _update_defaults(self, which):
        raise NotImplementedError


class AdvancedSubSettings(ipw.VBox):
    identifier = "sub"

    def __init__(self, model: AdvancedSubModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        self.loading_message = LoadingWidget(f"Loading {self.identifier} settings")

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_override_change,
            "override",
        )

        self.links = []

        self.rendered = False
        self.updated = False

    def render(self):
        raise NotImplementedError

    def refresh(self, which):
        self.updated = False
        self._unsubscribe()
        self._update(which)
        if hasattr(self._model, "input_structure") and not self._model.input_structure:
            self._model.reset()

    def _on_override_change(self, change):
        if not change["new"]:
            self._model.reset()

    def _unsubscribe(self):
        for link in self.links:
            link.unlink()
        self.links.clear()

    def _update(self, which):
        if self.updated:
            return
        self._model.update(which)
        self.updated = True
