import base64
from tempfile import NamedTemporaryFile
from threading import Thread, Event, Lock
from queue import Queue, Empty
from time import sleep
from collections import defaultdict
from itertools import count

import traitlets
import ipywidgets as ipw

from aiida.orm import CalcJobNode, ProcessNode


__all__ = [
    'ProcessOutputFollower',
]


_EOF = None


def _get_calcjobs(process):
    for called in process.called:
        if isinstance(called, CalcJobNode):
            yield called

        elif isinstance(called, ProcessNode):
            yield from _get_calcjobs(called)

        else:
            raise TypeError(called)


def _get_output(calcjob):
    assert isinstance(calcjob, CalcJobNode)
    if 'remote_folder' in calcjob.outputs:
        try:
            fn_out = calcjob.attributes['output_filename']
            with NamedTemporaryFile() as tmpfile:
                calcjob.outputs.remote_folder.getfile(fn_out, tmpfile.name)
                return tmpfile.read().decode().splitlines()
        except OSError as error:
            return list()
    else:
        return list()


class LogOutputWidget(ipw.VBox):

    value = traitlets.Unicode()
    value = traitlets.Tuple(traitlets.Unicode(), traitlets.Unicode())

    def __init__(self, num_min_lines=10, **kwargs):
        self._num_min_lines = num_min_lines

        self._filename = ipw.Text(
            description='Filename:',
            placeholder='No file selected.',
            disabled=True,
            layout=ipw.Layout(flex='1 1 auto', width='auto'),
        )

        self._output = ipw.HTML(layout=ipw.Layout(min_width='60em'))
        self._output_container = ipw.VBox(children=[self._output])
        self._download_link = ipw.HTML(layout=ipw.Layout(width='auto'))

        self._refresh_output()
        super().__init__(
            children=[
                self._output_container,
                ipw.HBox(children=[self._filename, self._download_link], layout=ipw.Layout(min_height='25px')),
            ], **kwargs)

    @traitlets.default('value')
    def _default_value(self):
        if self._num_min_lines > 0:
            return '', '\n' * self._num_min_lines

    @traitlets.observe('value')
    def _refresh_output(self, _=None):
        filename, loglines = self.value
        with self.hold_trait_notifications():
            self._filename.value = filename
            self._output.value = self._format_output(loglines)

            payload = base64.b64encode(loglines.encode()).decode()
            html_download = f'<a download="{filename}" href="data:text/plain;base64,{payload}" target="_blank">Download</a>'
            self._download_link.value = html_download

    def _format_output(self, text):
        lines = text.splitlines()

        # Add empty lines to reach the minimum number of lines.
        lines += [''] * max(0, self._num_min_lines - len(lines))

        # Replace empty lines with single white space to ensure that they are actually shown.
        lines = [line if len(line) > 0 else ' ' for line in lines]

        # Replace the first line if there is no output whatsoever
        if len(text.strip()) == 0 and len(lines) > 0:
            lines[0] = '[no output]'

        text = '\n'.join(lines)
        return '<pre style="background-color: #1f1f2e; color: white; line-height: 100%">{}</pre>'.format(text)


class ProcessOutputFollower(ipw.VBox):

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        # Setup the output queue and widgets.
        self._calcjobs = dict()

        self._output_queue = Queue()
        self._output = defaultdict(list)
        self._most_recently_updated_calcjob = None

        self._lock = Lock()
        self._push_threads = dict()
        self._stop_monitor_process = Event()
        self._pull_thread = None
        self._monitor_process_thread = None

        # This is the widget showing the various outputs.
        self._log_output = LogOutputWidget()

        # Setup selector which the user can use to select the displayed output node.
        self.selector = ipw.Dropdown(description='Calcjob:', layout=ipw.Layout(width='auto', flex='1 1 auto'))
        # Update the output widget when a different calcjob is selected.
        self.selector.observe(self._update, names=['value'])

        super().__init__(children=[self.selector, self._log_output], **kwargs)

    @traitlets.observe('process')
    def _observe_process(self, change):
        # Process any remaining items on the queue until it is empty.
        try:
            if change['old'].pk != change['new'].pk:
                # Old and new process are identical.
                return
        except AttributeError:
            pass

        with self._lock:
            self._stop_monitor_process.set()

            if self._monitor_process_thread:
                self._monitor_process_thread.join()
                self._monitor_process_thread = None

            if change['new']:
                # Initial update of all widgets.
                self._stop_monitor_process.clear()
                self._update()

                self._monitor_process_thread = Thread(target=self._monitor_process, args=(change['new'], ))
                self._monitor_process_thread.start()

    def _update_calcjobs(self, process):
        """Update the list of available calcjobs."""
        self._calcjobs = {cj.pk: cj for cj in _get_calcjobs(process)}

        with self.hold_trait_notifications():
            previous_calcjobs = {pk for _, pk in self.selector.options if pk is not None}
            new_calcjobs = set(self._calcjobs)

            if new_calcjobs != previous_calcjobs:
                options = [('<latest>', self._SHOW_MOST_RECENTLY_UPDATED)] + \
                    [(str(self._calcjobs[pk]), pk) for pk in sorted(self._calcjobs)]

                previously_selected = self.selector.value
                self.selector.options = options
                if previously_selected in self.selector.options:
                    self.selector.value = previously_selected

    def _start_push_threads(self):
        """Start push threads for calcjobs where thread is not started yet."""
        for pk, calcjob in self._calcjobs.items():
            if pk not in self._push_threads:
                push_thread = Thread(target=self._push_output, args=(calcjob, ))
                self._push_threads[pk] = push_thread
                push_thread.start()

    def _monitor_process(self, process):
        """Monitor process and orchestrate pushing and pulling of output."""
        i = count()
        while not self._stop_monitor_process.is_set():
            if next(i) % 10 == 0:
                self._update_calcjobs(process)
                self._start_push_threads()

            if self.process.is_sealed:
                break
            else:
                # Pull all new output and update widget.
                self._pull_output(timeout=.1)
                self._update()

        # One final search, push, pull, and update.
        self._update_calcjobs(process)
        self._start_push_threads()
        self._pull_output(timeout=5)
        self._update()

    def _push_output(self, calcjob, delay=.2):
        """Push new log lines onto the queue."""
        lineno = 0
        while True:
            try:
                lines = _get_output(calcjob)
            except Exception as error:
                self._output_queue.put((calcjob, [f'ERROR: {error}']))
            else:
                self._output_queue.put((calcjob, lines[lineno:]))
                lineno = len(lines)
            finally:
                if self._stop_monitor_process.is_set() or calcjob.is_sealed:
                    self._output_queue.put((calcjob, _EOF))
                    break
            sleep(delay)

    def _pull_output(self, timeout=None):
        """Pull new log lines from the queue and update output widget."""
        while True:
            try:
                calcjob, item = self._output_queue.get(timeout=timeout)
            except Empty:
                break
            else:
                if item is _EOF:
                    self._push_threads[calcjob.pk].join()
                else:
                    self._output[calcjob.pk].extend(item)
                self._most_recently_updated_calcjob = calcjob.pk

    _SHOW_MOST_RECENTLY_UPDATED = -1
    """Internal constant to signal that the widget should show the most recently updated output."""

    def _update(self, change=None):
        """Update the output widgets."""
        if change and (change['old'] == change['new']):
            return

        with self.hold_trait_notifications():
            pk = self.selector.value if change is None else change['new']

            if pk is self._SHOW_MOST_RECENTLY_UPDATED:
                calcjob = self._calcjobs.get(self._most_recently_updated_calcjob)
                restrict_num_lines = -10
            else:
                calcjob = self._calcjobs.get(pk)
                restrict_num_lines = None

            if calcjob is not None:
                filename = calcjob.attributes['output_filename']
                lines = self._output[calcjob.pk][restrict_num_lines:]
                self._log_output.value = filename, '\n'.join(lines)
