#!/usr/bin/env python3

"""This helper script is for temporarily remove the
`aiidalab-qe-workchain` package from dependencies list so that
in the test it will not use the remote released source.
"""

from configparser import ConfigParser

cf = ConfigParser()
cf.read("setup.cfg")
lst = cf["options"]["install_requires"].split("\n")
requirements_lst = [i for i in lst if not i.startswith("aiidalab-qe-workchain")]

with open("/tmp/requirements.txt", "w") as fh:
    fh.write("\n".join(requirements_lst[1:]))
