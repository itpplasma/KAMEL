"""Plot radial-current profiles for symmetric cyclotron-harmonic cutoffs."""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_PATH = Path(__file__).with_name("radial-current-harmonic-profiles.csv")
OUTPUT_STEM = Path(__file__).with_name("radial-current-harmonic-profiles")
MAX_CUTOFF = 6
RESONANT_RADIUS_CM = 36.8


def load_profiles(path: Path) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    expected_columns = {"radius_cm"}
    for cutoff in range(MAX_CUTOFF + 1):
        expected_columns.update({f"jrad_re_N{cutoff}", f"jrad_im_N{cutoff}"})

    missing_columns = expected_columns.difference(data.dtype.names or ())
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"missing profile columns: {missing}")

    radii = np.asarray(data["radius_cm"])
    profiles = {
        cutoff: np.asarray(data[f"jrad_re_N{cutoff}"]) + 1j * np.asarray(data[f"jrad_im_N{cutoff}"])
        for cutoff in range(MAX_CUTOFF + 1)
    }
    if not np.all(np.isfinite(radii)) or not all(
        np.all(np.isfinite(profile)) for profile in profiles.values()
    ):
        raise ValueError("profile data contain non-finite values")
    if not np.all(np.diff(radii) > 0.0):
        raise ValueError("profile radii must be strictly increasing")
    return radii, profiles


def main() -> None:
    radii, profiles = load_profiles(DATA_PATH)
    colors = plt.colormaps["viridis"](np.linspace(0.05, 0.95, MAX_CUTOFF + 1))

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8.5,
            "figure.dpi": 160,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), constrained_layout=True)

    ax = axes[0]
    for cutoff, color in zip(range(MAX_CUTOFF + 1), colors, strict=True):
        linewidth = 2.1 if cutoff in (0, MAX_CUTOFF) else 1.3
        linestyle = "--" if cutoff == 0 else "-"
        ax.plot(
            radii,
            np.abs(profiles[cutoff]),
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=rf"$N={cutoff}$",
        )
    ax.set_title(r"Cumulative profiles, $-N\leq\ell\leq N$")
    ax.set_ylabel(r"$|j^r_N|$ (statA cm$^{-2}$)")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax.legend(frameon=False, ncol=2)

    ax = axes[1]
    for cutoff, color in zip(range(1, MAX_CUTOFF + 1), colors[1:], strict=True):
        shell = profiles[cutoff] - profiles[cutoff - 1]
        ax.semilogy(
            radii,
            np.abs(shell),
            color=color,
            linewidth=1.7,
            label=rf"$\ell=\pm{cutoff}$",
        )
    ax.set_title("Added symmetric harmonic shells")
    ax.set_ylabel(r"$|j^r_N-j^r_{N-1}|$ (statA cm$^{-2}$)")
    ax.legend(frameon=False, ncol=2)

    for ax in axes:
        ax.set_xlabel(r"radius $r$ (cm)")
        ax.axvline(RESONANT_RADIUS_CM, color="0.35", linestyle=":", linewidth=1.0)
        ax.grid(which="major", color="0.84", linewidth=0.7)
        ax.grid(which="minor", color="0.92", linewidth=0.45)

    axes[0].annotate(
        r"$r_m=36.8$ cm",
        xy=(RESONANT_RADIUS_CM, 0.97),
        xycoords=("data", "axes fraction"),
        xytext=(4, 0),
        textcoords="offset points",
        ha="left",
        va="top",
        color="0.25",
        fontsize=8.5,
    )
    fig.suptitle(
        r"KIM radial-current profiles: static D, $(m,n)=(-6,2)$",
        fontsize=12,
    )
    for extension in ("png", "pdf"):
        fig.savefig(OUTPUT_STEM.with_suffix(f".{extension}"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
