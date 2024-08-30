# Stylesheets for the Quantum ESPRESSO app

This folder contains `css` stylesheets which may be loaded directly from the CSS folder using

```python
from aiidalab_widgets_base.utils.loaders import load_css_stylesheet

load_css_stylesheet(package="aiidalab_qe.app.static.styles.css", filename="<stylesheet-name>.css")
```

If `filename` is not provided, all CSS stylesheets will be loaded from the `css` package.
