"""Plot the radial-current cyclotron-harmonic convergence study."""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_DIR = Path(__file__).resolve().parent
CUTOFF = np.arange(1, 7)

# ||K_N - K_(N-1)||_F / ||K_6||_F
MATRIX_SHELLS = {
    r"$K^{j^r\Phi}$": np.array(
        [2.162318e-2, 1.641495e-3, 1.019670e-4, 4.952310e-6, 1.943242e-7, 6.373605e-9]
    ),
    r"$K^{j^r B^r}$": np.array(
        [3.860048e-9, 7.844008e-11, 2.360592e-12, 7.111152e-14, 1.956488e-15, 1.098316e-16]
    ),
    r"$K^{j^r B_\parallel}$": np.array(
        [1.112077e-3, 4.699163e-4, 3.247613e-5, 2.237533e-6, 1.140000e-7, 4.385394e-9]
    ),
}

# ||j_N - j_(N-1)||_2 / ||j_6||_2 and ||j_6 - j_N||_2 / ||j_6||_2
JRAD_SHELL = np.array(
    [2.274028e-2, 1.446188e-3, 9.229654e-5, 4.558101e-6, 1.804682e-7, 5.998154e-9]
)
JRAD_TAIL = np.array([1.541936e-3, 9.703313e-5, 4.744548e-6, 1.864662e-7, 5.998154e-9, np.nan])


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
            "figure.dpi": 160,
        }
    )

    colors = ["#0072B2", "#D55E00", "#009E73"]
    markers = ["o", "s", "^"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)

    ax = axes[0]
    for (label, values), color, marker in zip(MATRIX_SHELLS.items(), colors, markers, strict=True):
        ax.semilogy(
            CUTOFF,
            100.0 * values,
            marker=marker,
            color=color,
            linewidth=1.8,
            markersize=5,
            label=label,
        )
    ax.set_title("Kernel-matrix shell contribution")
    ax.set_ylabel("relative change (%)")
    ax.legend(frameon=False)

    ax = axes[1]
    ax.semilogy(
        CUTOFF,
        100.0 * JRAD_SHELL,
        "o-",
        color="#7B3294",
        linewidth=1.8,
        markersize=5,
        label=r"shell $\|j^r_N-j^r_{N-1}\|_2$",
    )
    ax.semilogy(
        CUTOFF,
        100.0 * JRAD_TAIL,
        "D--",
        color="#008837",
        linewidth=1.6,
        markersize=4.5,
        label=r"tail $\|j^r_6-j^r_N\|_2$",
    )
    ax.set_title("Reconstructed radial-current convergence")
    ax.legend(frameon=False)

    for ax in axes:
        ax.set_xlabel(r"symmetric cutoff $N$ ($-N\leq\ell\leq N$)")
        ax.set_xticks(CUTOFF)
        ax.grid(which="major", color="0.84", linewidth=0.7)
        ax.grid(which="minor", color="0.92", linewidth=0.45)
        ax.axhline(1.0, color="0.4", linewidth=0.8, linestyle=":")
        ax.axhline(0.1, color="0.55", linewidth=0.8, linestyle=":")

    fig.suptitle(
        r"KIM radial-current cyclotron harmonics: static D, $(m,n)=(-6,2)$",
        fontsize=12,
    )
    for extension in ("png", "pdf"):
        fig.savefig(
            OUTPUT_DIR / f"radial-current-harmonic-convergence.{extension}",
            bbox_inches="tight",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
