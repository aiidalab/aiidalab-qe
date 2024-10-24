import ipywidgets as ipw

from .model import AdvancedModel


class AdvancedSubSettings(ipw.VBox):
    identifier = "sub"

    def __init__(self, model: AdvancedModel, **kwargs):
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

    def _on_override_change(self, change):
        if not change["new"]:
            getattr(self._model, self.identifier).reset()

    def _unsubscribe(self):
        for link in self.links:
            link.unlink()
        self.links.clear()

    def _refresh(self):
        self.updated = False
        self._unsubscribe()
        self._update()

    def _update(self):
        raise NotImplementedError
