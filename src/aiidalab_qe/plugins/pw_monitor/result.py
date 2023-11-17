"""PW monitor result widget."""
import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultPanel


class Result(ResultPanel):
    title = "PW monitor"
    workchain_labels = []

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)

    def _update_view(self):
        """Update the view of the widget."""
        import threading

        update_interval = 1

        self.update_view = threading.Thread(
            target=self._draw_scf_convergence, args=(update_interval,)
        )
        self.update_view.start()

    def _draw_scf_convergence(self, update_interval):
        """Draw the SCF convergence plot."""
        import time

        import plotly.graph_objects as go

        #
        def get_data(relax_node):
            energies = []
            pw_base_nodes = [
                link.node
                for link in relax_node.base.links.get_outgoing().all()
                if link.link_label.startswith("iteration")
            ]
            for pw_base_node in pw_base_nodes:
                for link in pw_base_node.base.links.get_outgoing().all():
                    if not link.link_label.startswith("iteration"):
                        continue
                    energies.extend(link.node.base.extras.get("monitor_energies", []))
            return energies

        #
        if "relax" not in self.node.base.links.get_outgoing().all_link_labels():
            self.children = ipw.HTML("<h4>There is no relax calculation!</h4>")
            return
        # init figure
        relax_node = self.node.base.links.get_outgoing().get_node_by_label("relax")
        g = go.FigureWidget(
            layout=go.Layout(
                title=dict(text="SCF convergence"),
                barmode="overlay",
            )
        )
        g.layout.xaxis.title = "SCF step"
        g.layout.yaxis.title = "Energy (eV)"
        energies = get_data(relax_node)
        g.add_scatter(x=list(range(len(energies))), y=energies, mode="lines+markers")
        self.children = (ipw.HTML("<h4>PW calculation running!</h4>"), g)
        while relax_node.process_state.value not in ["finished", "excepted", "killed"]:
            energies = get_data(relax_node)
            with g.batch_update():
                g.data[0].x = list(range(len(energies)))
                g.data[0].y = energies
            time.sleep(update_interval)
        self.children = (ipw.HTML("<h4>PW calculation completed!</h4>"), g)

    def reset(self):
        """Reset the widget.
        Currently, this is not possible."""
        self.children = []
        self.update_view.stop()
