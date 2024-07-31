# QL-Light
This code calculates the local bifurcation criterion defined as
$$
\left.\frac{D^{\rm{ql}}_{e,22}}{D^{\rm{a}}}\right|_{r=r_s} \geq 1,
$$
where $D^{\rm{ql}}_{e,22}$ is the quasilinear heat diffusion coefficient of electrons [Markl NF 2023], $D^{\rm{a}}$ is the anomalous diffusion coefficient and $r_s$ is the radius of the resonant surface.
The threshold value on the right hand side is a heuristic value.

The quasilinear heat diffusion coefficient is regularly calculated in QL-Balance. However, to reduce the numerical effort and make the computation more practicable, QL-Light only calculates the necessary contributions to the diffusion coefficient.

