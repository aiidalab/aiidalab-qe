# -*- coding: utf8 -*-
from setuptools import setup, find_packages
import json

if __name__ == '__main__':

    with open('setup.json', 'r') as info:
        kwargs = json.load(info)

#    with open('requirements.txt', 'r') as rfile:
#        requirements = rfile.read().splitlines()

    setup(
        packages=find_packages(),
        include_package_data=True,
        reentry_register=True,
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        data_files=[
            # like `jupyter nbextension enable --sys-prefix`
            ("etc/jupyter/nbconfig/notebook.d", [
                "jupyter-config/nbconfig/notebook.d/aiidalab.json"
            ]),
            # like `jupyter serverextension enable --sys-prefix`
            ("etc/jupyter/jupyter_notebook_config.d", [
                "jupyter-config/jupyter_notebook_config.d/aiidalab.json"
            ])
        ],
#        install_requires=requirements,
        # For packages not yet published on pypi
        # Note: this only works when installing aiidalab from github
        zip_safe=False,
        **kwargs
    )
