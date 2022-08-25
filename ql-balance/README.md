# ql-balance
To compile the code, run

	make -f Balance.mk_mpi

This will compile the code and create the executable balance.x.mpif90.openmpi_x86_64.

The file balance_conf.nml is the blueprint namelist file required by the ql-balance code.

# Run balance code
Use

LD_LIBRARY_PATH=/proj/plasma/soft/math_libs/64bit/sundials-2.6.2/lib/ ./balance.x.mpif90.openmpi_x86_64

in the run directory.
