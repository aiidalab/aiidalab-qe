"""Widgets for the monitoring of processes."""
import os
from collections import deque

import traitlets
import ipywidgets as ipw
from aiida.orm import ProcessNode


def get_calc_job_output(process):
    from aiidalab_widgets_base.process import get_running_calcs
    previous_calc_id = None
    num_lines = 0

    while not process.is_sealed:
        calc = None
        for calc in get_running_calcs(process):
            if calc.id == previous_calc_id:
                break
        else:
            if calc:
                previous_calc_id = calc.id

        if calc and 'remote_folder' in calc.outputs:
            f_path = os.path.join(calc.outputs.remote_folder.get_remote_path(),
                                  calc.attributes['output_filename'])
            if os.path.exists(f_path):
                with open(f_path) as fobj:
                    new_lines = fobj.readlines()[num_lines:]
                    num_lines += len(new_lines)
                    yield from new_lines


class ProgressBarWidget(ipw.VBox):
    """A bar showing the proggress of a process."""

    def __init__(self, **kwargs):
        self._process = None
        self.correspondance = {
            None: (0, 'warning'),
            "created": (0, 'info'),
            "running": (1, 'info'),
            "waiting": (1, 'info'),
            "killed": (2, 'danger'),
            "excepted": (2, 'danger'),
            "finished": (2, 'success'),
        }
        self.bar = ipw.IntProgress(  # pylint: disable=blacklisted-name
            value=0,
            min=0,
            max=2,
            step=1,
            bar_style='warning',  # 'success', 'info', 'warning', 'danger' or ''
            orientation='horizontal',
            layout=ipw.Layout(width="auto")
        )
        self.state = ipw.HTML(
            description="Calculation state:", value='',
            style={'description_width': '100px'},
        )
        super().__init__(children=[self.state, self.bar], **kwargs)

    def update(self):
        """Update the bar."""
        self.bar.value, self.bar.bar_style = self.correspondance[self.current_state]
        if self.current_state is None:
            self.state.value = 'N/A'
        else:
            self.state.value = self.current_state.capitalize()

    @property
    def process(self):
        return self._process

    @process.setter
    def process(self, value):
        self._process = value
        self.update()

    @property
    def current_state(self):
        if self.process is not None:
            return self.process.process_state.value


class LogOutputWidget(ipw.VBox):

    def __init__(self, title='Output:', num_lines_shown=3, **kwargs):
        self.description = ipw.Label(value=title)
        self.last_lines = ipw.HTML()

        self.lines = []
        self.lines_shown = deque([''] * num_lines_shown, maxlen=num_lines_shown)

        self.raw_log = ipw.Textarea(
            layout=ipw.Layout(
                width='auto',
                height='auto',
                display='flex',
                flex='1 1 auto',
            ),
            disabled=True)

        self.accordion = ipw.Accordion(children=[self.raw_log])
        self.accordion.set_title(0, 'Raw log')
        self.accordion.selected_index = None

        self._update()
        super().__init__(children=[self.description, self.last_lines, self.accordion], **kwargs)

    def clear(self):
        self.lines.clear()
        self.lines_shown.clear()
        self._update()

    def append_line(self, line):
        self.lines.append(line.strip())
        self.lines_shown.append("{:03d}: {}".format(len(self.lines), line.strip()))
        self._update()

    @staticmethod
    def _format_code(text):
        return '<pre style="background-color: #1f1f2e; color: white;">{}</pre>'.format(text)

    def _update(self):
        with self.hold_trait_notifications():
            lines_to_show = self.lines_shown.copy()
            while len(lines_to_show) < self.lines_shown.maxlen:
                lines_to_show.append(' ')

            self.last_lines.value = self._format_code('\n'.join(lines_to_show))
            self.raw_log.value = '\n'.join(self.lines)


class ProcessStatusWidget(ipw.VBox):

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        self._process_output = None  # Generator function for the process output

        # Widgets
        self.progress_bar = ProgressBarWidget()
        self.log_output = LogOutputWidget()
        self.process_id_text = ipw.Text(
            value='',
            description='Process:',
            layout=ipw.Layout(width='auto', flex="1 1 auto"),
            disabled=True,
        )
        ipw.dlink((self, 'process'), (self.process_id_text, 'value'),
                  transform=lambda proc: str(proc))

        super().__init__(children=[self.progress_bar,
                                   self.log_output, self.process_id_text], **kwargs)

    @traitlets.observe('process')
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            self._process_output = get_calc_job_output(change['new'])
            self.progress_bar.process = change['new']
            if change['new'] != change['old']:
                self.log_output.clear()

    def update(self):
        try:
            self.log_output.append_line(next(self._process_output))
        except (TypeError, StopIteration):
            pass
        self.progress_bar.update()
