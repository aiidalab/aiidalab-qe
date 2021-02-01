import ipywidgets as ipw


def get_start_widget(appbase, jupbase, notebase):
    # http://fontawesome.io/icons/
    template = """
    <div align="center">
        <a href="{appbase}/qe.ipynb" target="_blank">
            <img src="https://gitlab.com/QEF/q-e/raw/develop/logo.jpg" height="120px" width=243px">
        </a>
    </div>
    """

    html = template.format(appbase=appbase, jupbase=jupbase, notebase=notebase)
    return ipw.HTML(html)


# EOF
