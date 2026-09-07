# Gauss segment-transform audit

`verify_gauss_segment_transforms.wl` is an executable symbolic derivation for
the standard-Gauss global-FEM `F1` and `F2` cell transforms. It complements
`test_gauss_segment_transforms.f90`, whose expected values come from direct
radial quadrature and numerical derivatives.

The derivation starts from the Fourier kernel in Chapter 13, especially the
definition of `b+` in Eq. (13.46), and applies the Fourier-to-hat-basis
transformation of Eqs. (15.13) and (15.17). With

```text
c = cos(theta), s = sin(theta)
mu = (x + x') / 2, delta = x - x'
```

the radial Gaussian generator is

```text
Q = exp[-ks^2 rhoT^2
        -(rg-mu)^2/(rhoT^2 (1+c))
        -delta^2/(4 rhoT^2 (1-c))] / (rhoT^2 s).
```

The script proves the following identities for arbitrary finite cell bounds,
`rhoT > 0`, and `0 < theta < pi`:

1. Directly integrating `Q` produces the `sqrt(1+cos(theta))` factor in `F1`.
2. The closed `Jrg1`, `Jrg23`, and `Jrg4` expressions are the corresponding
   Gaussian cell moments.
3. Applying the `b+` multiplier in real space,

   ```text
   Q2 = rhoT^2 ks^2 Q
        - rhoT^2/2 (d_x^2 + d_x'^2) Q,
   ```

   gives the corrected `F2` bracket

   ```text
   2 s^2 (1 + ks^2 rhoT^2 s^2) Jrg1
   - (1 + c^2) Jrg23 + 4 c Jrg4.
   ```

Run it with a licensed Wolfram kernel:

```bash
WolframKernel -script KIM/tests/symbolic/verify_gauss_segment_transforms.wl
```

The endpoint representations are singular and are interpreted by limits;
the production Gauss-Legendre rule does not evaluate the endpoints.
