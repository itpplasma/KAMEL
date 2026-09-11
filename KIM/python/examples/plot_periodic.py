"""Plot fields from an existing periodic KIM HDF5 result.

Usage::

    python examples/plot_periodic.py runs/<run-id>/output.h5 periodic.png

The script only reads the supplied result. It never invokes ``KIM.x``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from kim import Result


def plot_periodic_result(result_path: Path | str, figure_path: Path | str) -> Any:
    """Save a two-panel plot of a stored periodic result and return its closed figure.

    KIM writes the periodic potential as ``statV`` and parallel current density as
    ``statA/cm^2`` (``poisson_periodic.f90`` lines 627--634). The radial field grid
    is the effective radius in centimetres. The returned figure is already closed
    so callers can inspect it without leaving an interactive figure open.
    """

    # Keep matplotlib optional: importing the KIM API and running the CLI does not
    # require the plotting extra.
    from matplotlib import pyplot as plt

    periodic = Result(result_path).periodic
    radius = periodic.radius
    resonance = periodic.resonance_radius
    half_width = periodic.as_is_half_width

    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(8.0, 6.0), constrained_layout=True)
    potential_axis, current_axis = axes
    potential = periodic.potential
    current = periodic.parallel_current_density

    potential_axis.plot(radius, potential.real, label="Re(Phi)")
    potential_axis.plot(radius, potential.imag, label="Im(Phi)")
    current_axis.plot(radius, current.real, label="Re(jpar)")
    current_axis.plot(radius, current.imag, label="Im(jpar)")

    for axis in axes:
        axis.axvline(resonance, color="black", linestyle="--", label="resonance")
        axis.axvspan(
            resonance - half_width,
            resonance + half_width,
            color="tab:orange",
            alpha=0.2,
            label="as-is interval",
        )
        axis.grid(True, alpha=0.25)
        axis.legend()

    potential_axis.set_ylabel("Phi [statV]")
    potential_axis.set_xlabel("r_eff [cm]")
    current_axis.set_ylabel("jpar [statA/cm^2]")
    current_axis.set_xlabel("r_eff [cm]")
    potential_axis.set_title("KIM forced-periodicity response")

    destination = Path(figure_path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        figure.savefig(destination)
    finally:
        plt.close(figure)
    return figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path, help="existing KIM HDF5 result")
    parser.add_argument("figure", type=Path, help="output image path")
    arguments = parser.parse_args()
    plot_periodic_result(arguments.result, arguments.figure)


if __name__ == "__main__":
    main()
