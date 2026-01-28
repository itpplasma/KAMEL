import os

import f90nml
import numpy as np


class balance_conf:

    blueprint_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "QL-Balance", "namelists", "balance_conf.nml"
    )

    def __init__(self, path=""):
        if path == "":
            self.conf = f90nml.read(self.blueprint_path)
        else:
            self.conf = f90nml.read(path)

    def write_conf(self, path):
        self.conf.write(path, force=True)
