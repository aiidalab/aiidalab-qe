# Stylesheets for the Quantum ESPRESSO app

This folder contains two folders:

- `scss` - a folder for SCSS-format stylesheets
- `css`  - a folder serving as the compilation target for SCSS files

It is recommended to add SCSS stylesheets in the SCSS folder for functionality and flexibility (see https://sass-lang.com/).

Compiling SCSS into CSS (in the CSS folder) can be done using:

- The [compile-hero](https://marketplace.visualstudio.com/items?itemName=Wscats.eno) extension, if using VS Code
- The `sass` [CLI](https://sass-lang.com/install/)

Stylesheets may be loaded directly from the CSS folder using

```python
from aiidalab_widgets_base.utils.loaders import load_css_stylesheet

load_css_stylesheet(package="aiidalab_qe.app.static.styles.css", filename="<stylesheet-name>.css")
```

If `filename` is not provided, all CSS stylesheets will be loaded from the `css` package.
