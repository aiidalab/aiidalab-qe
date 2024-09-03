# Stylesheets for the Quantum ESPRESSO app

This folder contains `css` stylesheets. These can be loaded from the styles folder using

```python
from aiidalab_widgets_base.utils.loaders import load_css

load_css(css_path="src/aiidalab_qe/app/static/styles")  # load all stylesheets in the styles folder

# or

load_css(css_path="src/aiidalab_qe/app/static/styles/<stylesheet>.css") # load a single stylesheet
```
