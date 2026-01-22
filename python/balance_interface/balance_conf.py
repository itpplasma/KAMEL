import numpy as np
import os

try:
    import f90nml
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "f90nml is required for balance_interface. Install it with "
        "`python3 -m pip install f90nml`."
    ) from exc

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
