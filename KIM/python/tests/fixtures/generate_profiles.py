"""Deterministic parabolic `r_eff` profiles used by Python API tests."""

from pathlib import Path

import numpy as np


def generate_profiles(directory: Path) -> dict[str, np.ndarray]:
    """Generate the reference profiles; physical values use KIM's CGS conventions."""

    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)

    radius = np.linspace(0.0, 10.0, 11)
    normalized_radius = radius / radius[-1]
    parabolic_shape = 1.0 - normalized_radius**2
    profiles = {
        "radius": radius,
        "n": 1.0e13 + (4.5e13 - 1.0e13) * parabolic_shape,  # 1/cm^3
        "Te": 600.0 + (4_000.0 - 600.0) * parabolic_shape,  # eV
        "Ti": 600.0 + (3_000.0 - 600.0) * parabolic_shape,  # eV
        # q=-3.5 at r_eff=5 cm: the (m,n)=(7,2) resonant surface.
        "q": -(3.5 + 0.1 * (radius - 5.0)),
        # Signed quantities intentionally remain negative through the domain.
        "Er": -0.1 - 0.4 * parabolic_shape,  # statV/cm
        "Vz": -50_000.0 - 200_000.0 * parabolic_shape,  # cm/s
    }
    for name in ("n", "Te", "Ti", "q", "Er", "Vz"):
        np.savetxt(
            destination / f"{name}.dat",
            np.column_stack((radius, profiles[name])),
            fmt="%.16e",
        )
    return profiles


if __name__ == "__main__":
    generate_profiles(Path(__file__).with_name("profiles"))
