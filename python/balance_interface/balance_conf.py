import math
import os
from numbers import Integral, Real

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
        keeps the policy explicit: drift-kinetic electrons,
        finite-Larmor-radius ions, and no compression drive until the coupled
        B-parallel response is implemented.
        """
        modes = list(modes)
        if not modes:
            raise ValueError("at least one nonzero (m, n) mode is required")
        if len(modes) > 100:
            raise ValueError("periodic KIM supports at most 100 modes")
        normalized_modes = []
        for mode in modes:
            if isinstance(mode, (str, bytes)):
                raise ValueError("periodic KIM modes must be nonzero (m, n) pairs")
            try:
                m_mode, n_mode = mode
            except (TypeError, ValueError) as exc:
                raise ValueError("periodic KIM modes must be nonzero (m, n) pairs") from exc
            if (
                any(
                    isinstance(value, bool) or not isinstance(value, Integral)
                    for value in (m_mode, n_mode)
                )
                or m_mode == 0
                or n_mode == 0
            ):
                raise ValueError("periodic KIM modes must be nonzero (m, n) pairs")
            normalized_modes.append((int(m_mode), int(n_mode)))
        if len(set(normalized_modes)) != len(normalized_modes):
            raise ValueError("periodic KIM modes must be unique")
        if (
            isinstance(target_current, bool)
            or not isinstance(target_current, Real)
            or not math.isfinite(target_current)
            or target_current < 0.0
        ):
            raise ValueError("target_current must be a finite non-negative number")
        if benchmark_mode not in ("none", "drift_kinetic_limit"):
            raise ValueError("benchmark_mode must be none or drift_kinetic_limit")

        balance = self.conf.setdefault("balancenml", {})
        balance.update(
            {
                "wave_code": "KIM",
                "kim_run_type": "electrostatic_periodic",
                "kim_profiles_from_balance": True,
                "kim_electron_transport_model": "drift_kinetic",
                "kim_ion_transport_model": "finite_larmor_radius",
                "kim_bparallel_source": "disabled",
                "kim_benchmark_mode": benchmark_mode,
                "kim_transport_benchmark": benchmark_mode == "drift_kinetic_limit",
                "kim_n_modes": len(modes),
                "kim_m_list": [mode[0] for mode in normalized_modes],
                "kim_n_list": [mode[1] for mode in normalized_modes],
                "I_par_toroidal": float(target_current),
            }
        )
