# Radial-current cyclotron-harmonic convergence

## Question

Do nonzero cyclotron harmonics materially change KIM's radial-current row,
including the reconstructed radial current density?

## Case and method

The convergence harness `test_radial_current_harmonics` uses a static
Fokker--Planck deuterium plasma with (B_0=-18{,}000\) G and mode
((m,n)=(-6,2)). The analytic safety-factor profile crosses (q=3) at
(r_m=36.8\) cm. The periodic window uses 64 gyrocenter points, Fourier cutoff
(M=3), and as-is/transition half-widths of 5/10 ion Larmor radii.

For a symmetric harmonic cutoff (N), each kernel sums
(ell=-N,\ldots,+N). Matrix-shell changes are measured as

\[
 \frac{\lVert K_N-K_{N-1}\rVert_F}{\lVert K_6\rVert_F}.
\]

The reconstructed-current shell and remaining tail use the discrete radial
(L_2) norm, normalized by (lVert j^r_6\rVert_2). The potential is solved
once; for each cutoff,

\[
 j^r_N=K^{j^r\Phi}_N\Phi+K^{j^rB^r}_N B^r.
\]

The implemented (K^{j^rB_\parallel}) is measured independently but is not
included in (j^r_N), because the present periodic solver has no
(B_\parallel) field.

## Results

| Shell | (K^{j^r\Phi}) | (K^{j^rB^r}) | (K^{j^rB_\parallel}) | (j^r) shell | (j^r) tail after shell |
|---:|---:|---:|---:|---:|---:|
| 1 | 2.162318e-2 | 3.860048e-9 | 1.112077e-3 | 2.274028e-2 | 1.541936e-3 |
| 2 | 1.641495e-3 | 7.844008e-11 | 4.699163e-4 | 1.446188e-3 | 9.703313e-5 |
| 3 | 1.019670e-4 | 2.360592e-12 | 3.247613e-5 | 9.229654e-5 | 4.744548e-6 |
| 4 | 4.952310e-6 | 7.111152e-14 | 2.237533e-6 | 4.558101e-6 | 1.864662e-7 |
| 5 | 1.943242e-7 | 1.956488e-15 | 1.140000e-7 | 1.804682e-7 | 5.998154e-9 |
| 6 | 6.373605e-9 | 1.098316e-16 | 4.385394e-9 | 5.998154e-9 | 0 (reference) |

The (N=0\) radial current differs from the (N=6\) reference by
2.394018e-2. The reference current norm is
(lVert j^r_6\rVert_2=8.741477\times10^8\) in the test's discrete CGS
normalization.

## Interpretation

The nonzero harmonics are relevant for (j^r) in this case. Retaining only
(ell=0) misses about 2.4% of the reconstructed current. Most of that change
comes from the (ell=\pm1) shell through (K^{j^r\Phi}); the
(ell=\pm2) shell remains a 0.14% correction. A cutoff of (N=3) leaves a
4.7e-6 relative tail, while (N=4) leaves 1.9e-7, so (N=4) is ample for
this case.

The three columns behave differently. (K^{j^rB^r}) is insensitive to
nonzero harmonics here, whereas (K^{j^r\Phi}) has a percent-level first
shell. (K^{j^rB_\parallel}) has the expected nonzero neighboring-harmonic
structure, but its quantitative effect on total (j^r) cannot be assessed
until a (B_\parallel) field is supplied.

These numbers characterize this static deuterium case, not a universal
ordering. In particular, finite-frequency cases near a cyclotron resonance
require a fresh cutoff scan.
