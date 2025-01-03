# QL-Balance
QL-Balance calculates the quasilinear plasma response, i.e. the impact of magnetic perturbations on the plasma profiles, within a set of transport equations. The quasilinear transport coefficients are computed with the EM perturbation fields provided by KiLCA. Relevant publications are Heyn NF 2014 and Markl NF 2023.

To compile the code, use
```
make
```

This will compile the code and create the executables ql-balance and ql-light in the build directory.

The file balance_conf.nml is the blueprint namelist file required by the code.

Note that the prior build system was make with the make file Balance.mk_mpi. This created the executable balance.x.mpif90.openmpi_x86_64.

# Run balance code
Running the balance code requires the balance_conf.nml file in the run directory, as well as the usual KiLCA input files.
