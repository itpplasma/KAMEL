OBJS =  OBJS/plag_coeff.o \
        OBJS/sample_matrix_mod.o \
        OBJS/get_matrix.o \
        OBJS/sample_matrix.o \
        OBJS/split_normal_modes.o \
        OBJS/lapack_interfs.o \
        OBJS/odeint_allroutines.o \
        OBJS/residual_system.o \
        OBJS/binsrc.o \
        OBJS/stiff_solver.o \
        OBJS/wave_stuff.o \
        OBJS/driver_separate_eigenmodes.o
COMP = f95-lah
OPTS= -O -M OBJS
#OPTS = --chk a,e,s,u,x -M OBJS
#COMP = g95
#OPTS =  -ftrace=full
#COMP =  gfortran-4.3.1
#OPTS =  -fbounds-check
#COMP = f95-intel
#OPTS =  -check all
separate_eigenmodes.x: $(OBJS) Separate_eigenmodes.mk
    $(COMP) $(OPTS) -o separate_eigenmodes.x $(OBJS) -llapack
OBJS/plag_coeff.o: plag_coeff.f90 Separate_eigenmodes.mk
    $(COMP) $(OPTS) -c plag_coeff.f90
    mv plag_coeff.o OBJS
OBJS/sample_matrix_mod.o: sample_matrix_mod.f90 Separate_eigenmodes.mk
    $(COMP) $(OPTS) -c sample_matrix_mod.f90
    mv sample_matrix_mod.o OBJS
OBJS/get_matrix.o: get_matrix.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c get_matrix.f90
    mv get_matrix.o OBJS
OBJS/sample_matrix.o: sample_matrix.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c sample_matrix.f90
    mv sample_matrix.o OBJS
OBJS/split_normal_modes.o: split_normal_modes.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c split_normal_modes.f90
    mv split_normal_modes.o OBJS
OBJS/lapack_interfs.o: lapack_interfs.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c lapack_interfs.f90
    mv lapack_interfs.o OBJS
OBJS/odeint_allroutines.o: odeint_allroutines.f Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c odeint_allroutines.f
    mv odeint_allroutines.o OBJS
OBJS/residual_system.o: residual_system.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c residual_system.f90
    mv residual_system.o OBJS
OBJS/binsrc.o: binsrc.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c binsrc.f90
    mv binsrc.o OBJS
OBJS/stiff_solver.o: stiff_solver.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c stiff_solver.f90
    mv stiff_solver.o OBJS
OBJS/wave_stuff.o: wave_stuff.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c wave_stuff.f90
    mv wave_stuff.o OBJS
OBJS/driver_separate_eigenmodes.o: driver_separate_eigenmodes.f90 Separate_eigenmodes.mk sample_matrix_mod.f90
    $(COMP) $(OPTS) -c driver_separate_eigenmodes.f90
    mv driver_separate_eigenmodes.o OBJS
#OBJS/.o: .f90 Separate_eigenmodes.mk sample_matrix_mod.f90
#   $(COMP) $(OPTS) -c .f90
#   mv .o OBJS
