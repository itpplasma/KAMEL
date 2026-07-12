#!/usr/bin/env python3
"""Compare the hat-basis and periodic-Fourier electrostatic potentials.

No-op benchmark for the enforced-periodicity solver (issue #175): with the
periodization layer covering the localized response region, the periodic
Fourier solution of the screened Poisson system must reproduce the hat-basis
solution there.  Both solves are driven by the same constant unit Br
(type_br_field = 12) on the same Krook background.

Inputs are the text outputs of the two runs:
  hat:      <hat_run>/out/m-3_n2/fields/Phi.dat           (r, Re, Im, |Phi|)
  periodic: <per_run>/out/m-3_n2/fields/Phi_periodic.dat  (r, Re, Im)

The comparison window is the flat part of the periodization layer,
|r - fp_r_res| <= fp_dr_layer, where the periodized background is identical
to the original one.  The hat solution is interpolated onto the periodic
output grid restricted to that window.  A k = 0 gauge constant can differ
between the two bases (the hat solve pins the boundary values, the periodic
solve fixes the mean over one period), so the metrics are reported both raw
and after removing the complex constant that minimizes the L2 difference.

Exit status 0 iff the offset-removed relative L2 difference is below the
threshold.  On failure the windowed profiles are written next to the outputs
as hat_window.txt and periodic_window.txt.
"""

import argparse
import sys

import numpy as np

THRESHOLD_REL_L2 = 0.05


def load_complex_profile(path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1] + 1j * data[:, 2]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("hat_phi", help="hat-basis fields/Phi.dat")
    parser.add_argument("periodic_phi",
                        help="periodic fields/Phi_periodic.dat")
    parser.add_argument("--r-res", type=float, default=19.995,
                        help="fp_r_res of the periodic run [cm]")
    parser.add_argument("--dr-layer", type=float, default=2.3,
                        help="fp_dr_layer of the periodic run [cm]")
    parser.add_argument("--threshold", type=float, default=THRESHOLD_REL_L2,
                        help="acceptance threshold on offset-removed "
                             "relative L2")
    parser.add_argument("--outdir", default=".",
                        help="where to write the windowed profiles on "
                             "failure")
    args = parser.parse_args()

    r_hat, phi_hat = load_complex_profile(args.hat_phi)
    r_per, phi_per = load_complex_profile(args.periodic_phi)

    window = np.abs(r_per - args.r_res) <= args.dr_layer
    r_cmp = r_per[window]
    per_cmp = phi_per[window]
    hat_cmp = (np.interp(r_cmp, r_hat, phi_hat.real)
               + 1j * np.interp(r_cmp, r_hat, phi_hat.imag))

    hat_norm_l2 = np.linalg.norm(hat_cmp)
    hat_norm_inf = np.abs(hat_cmp).max()

    def metrics(offset):
        diff = hat_cmp - per_cmp - offset
        return (np.linalg.norm(diff) / hat_norm_l2,
                np.abs(diff).max() / hat_norm_inf)

    raw_l2, raw_inf = metrics(0.0)
    gauge = np.mean(hat_cmp - per_cmp)
    off_l2, off_inf = metrics(gauge)

    print(f"comparison window: |r - {args.r_res}| <= {args.dr_layer} cm "
          f"({r_cmp.size} points, [{r_cmp[0]:.3f}, {r_cmp[-1]:.3f}] cm)")
    print(f"hat   max |Phi| in window: {np.abs(hat_cmp).max():.6e} statV")
    print(f"per   max |Phi| in window: {np.abs(per_cmp).max():.6e} statV")
    print(f"gauge constant (mean difference): {gauge:.6e}")
    print(f"raw            rel L2 = {raw_l2:.6e}   rel Linf = {raw_inf:.6e}")
    print(f"offset-removed rel L2 = {off_l2:.6e}   rel Linf = {off_inf:.6e}")

    if off_l2 <= args.threshold:
        print(f"PASS: offset-removed rel L2 {off_l2:.3e} "
              f"<= {args.threshold:.3e}")
        return 0

    header = "r [cm]  Re(Phi) [statV]  Im(Phi) [statV]"
    np.savetxt(f"{args.outdir}/hat_window.txt",
               np.column_stack([r_cmp, hat_cmp.real, hat_cmp.imag]),
               header=header)
    np.savetxt(f"{args.outdir}/periodic_window.txt",
               np.column_stack([r_cmp, per_cmp.real, per_cmp.imag]),
               header=header)
    print(f"FAIL: offset-removed rel L2 {off_l2:.3e} > {args.threshold:.3e}")
    print(f"windowed profiles written to {args.outdir}/hat_window.txt "
          f"and {args.outdir}/periodic_window.txt")
    return 1


if __name__ == "__main__":
    sys.exit(main())
