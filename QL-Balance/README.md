# QL-Balance
QL-Balance calculates the quasilinear plasma response, i.e. the impact of magnetic perturbations on the plasma profiles, within a set of transport equations. The quasilinear transport coefficients are computed with the EM perturbation fields provided by KiLCA. Relevant publications are Heyn NF 2014 and Markl NF 2023.

To compile the code, use
```
./ql-balance.sh
```

This will compile the code and create the executables ql-balance and ql-light in the build directory.

The file balance_conf.nml is the blueprint namelist file required by the code.

Note that the prior build system was make with the make file Balance.mk_mpi. This created the executable balance.x.mpif90.openmpi_x86_64.

# Run balance code
Use

LD_LIBRARY_PATH=/proj/plasma/soft/math_libs/64bit/sundials-2.6.2/lib/ ./balance.x.mpif90.openmpi_x86_64

in the run directory.


# QL-Light
Light version of the balance code that only calculates Dqle22 for the local bifurcation criterion. Also includes constant-psi approximation.

### todo:
- [ ] Read in data
- [ ] Determine gradients
- [ ] Calculate quasilinear diffusion coefficient Dqle22
- [ ] Calculate quasilinear diffusion coefficient Dqle22 in const-psi
- [ ] Optional: calculate torque
- [ ] Optional: calculate torque in const-psi
- [ ] interpolate results
- [ ] write out results