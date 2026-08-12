import os

import numpy as np

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

    def configure_periodic_kim(self, modes, target_current=0.0, benchmark_mode="none"):
        """Select the validated periodic-KIM transport contract.

        ``modes`` is an iterable of ``(m, n)`` pairs.  The generated namelist
        keeps the policy explicit: drift-kinetic electrons, integral ions, and
        the self-consistent periodic B-parallel source.
        """
        modes = list(modes)
        if not modes:
            raise ValueError("at least one nonzero (m, n) mode is required")
        if len(modes) > 100:
            raise ValueError("periodic KIM supports at most 100 modes")
        if any(len(mode) != 2 or int(mode[0]) == 0 or int(mode[1]) == 0 for mode in modes):
            raise ValueError("periodic KIM modes must be nonzero (m, n) pairs")
        if float(target_current) < 0.0:
            raise ValueError("target_current must be non-negative")
        if benchmark_mode not in ("none", "drift_kinetic_limit"):
            raise ValueError("benchmark_mode must be none or drift_kinetic_limit")

        balance = self.conf.setdefault("balancenml", {})
        balance.update(
            {
                "wave_code": "KIM",
                "kim_run_type": "electrostatic_periodic",
                "kim_profiles_from_balance": True,
                "kim_electron_transport_model": "drift_kinetic",
                "kim_ion_transport_model": "integral",
                "kim_bparallel_source": "periodic",
                "kim_benchmark_mode": benchmark_mode,
                "kim_n_modes": len(modes),
                "kim_m_list": [int(mode[0]) for mode in modes],
                "kim_n_list": [int(mode[1]) for mode in modes],
                "I_par_toroidal": float(target_current),
            }
        )
