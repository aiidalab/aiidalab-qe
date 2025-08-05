import traitlets as tl

from aiidalab_qe.common.mvc import Model

from ..utils import (
    FUNCTIONAL_MAPPING,
    UpfData,
    get_info_from_family,
    get_pseudo_family_by_content,
    get_upf_dict,
)


class PseudoPotentialUploaderModel(Model):
    pseudo = tl.Instance(UpfData, allow_none=True)
    filename = tl.Unicode(allow_none=True)
    info = tl.Unicode(allow_none=True)
    cutoffs = tl.List(tl.Float(), [])
    uploaded = tl.Bool(False)
    message = tl.Unicode(allow_none=True)

    PSEUDO_INFO_MESSAGE = """
        <div class="pseudo-text">
            {ecutwfc} | {ecutrho} | {functional} | {relativistic}
        </div>
    """

    def __init__(self, kind_name: str, kind_symbol: str, **kwargs):
        super().__init__(**kwargs)
        self.kind_name = kind_name
        self.kind_symbol = kind_symbol

    def update_pseudo_info(self):
        if not self.pseudo:
            self.info = ""
            return
        self.info = self.PSEUDO_INFO_MESSAGE.format(
            ecutwfc=self.cutoffs[0] or "?",
            ecutrho=self.cutoffs[1] or "?",
            functional=self.pseudo.base.extras.get("functional", None) or "?",
            relativistic=self.pseudo.base.extras.get("relativistic", None) or "N/A",
        )

    def populate_extras(self):
        if not self.pseudo:
            return

        keys = ("functional", "relativistic")
        if all(key in self.pseudo.base.extras.all for key in keys):
            return

        try:
            upf_dict = get_upf_dict(self.pseudo)
        except Exception as err:
            raise ValueError(
                f"Failed to read UPF data from pseudo potential with UUID {self.pseudo.uuid}"
            ) from err

        if pseudo_family := get_pseudo_family_by_content(self.pseudo.md5):
            functional, relativistic = get_info_from_family(pseudo_family)
        else:
            functional = FUNCTIONAL_MAPPING.get(
                functional := "+".join(
                    upf_dict["header"].get("functional", "").lower().split()
                ),
                functional,
            )
            relativistic = upf_dict["header"].get("relativistic", None)

        self.pseudo.base.extras.set_many(
            {
                "functional": functional or None,
                "relativistic": relativistic,
            }
        )
