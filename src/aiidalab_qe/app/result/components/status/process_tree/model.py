import traitlets as tl

from aiidalab_qe.common.mixins import HasProcess
from aiidalab_qe.common.mvc import Model


class SimplifiedProcessTreeModel(Model, HasProcess):
    clicked = tl.Unicode(None, allow_none=True)
