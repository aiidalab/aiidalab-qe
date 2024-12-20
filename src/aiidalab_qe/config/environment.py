from pathlib import Path

import yaml

CONFIG_PATH = Path.home() / ".aiidalab" / "config.yml"


def on_demo_server():
    if CONFIG_PATH.exists():
        with CONFIG_PATH.open("r") as file:
            config = yaml.safe_load(file)
            return config.get("on_demo_server", False)
