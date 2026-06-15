# Only include the actual DDEABM ODE-solver sources.
# The slatec error-handling helpers (i1mach, d1mach, j4save, xermsg, xerprn,
# xersve, xercnt, xerhlt, xgetua, fdump) are provided by libslatec and must
# NOT be duplicated here.  ddeabm's i1mach has all DATA statements commented
# out (returns 0 for every query), which breaks slatec's DGAMLN/Bessel
# functions when it shadows the correct slatec version at link time.
set(ddeabm_files
    ddeabm.f
    ddes.f
    dhstrt.f
    dhvnrm.f
    dintp.f
    dsteps.f
)
add_library(ddeabm "${ddeabm_files}")
# ddeabm.f/ddes.f reference xermsg (and i1mach/d1mach via that path); those
# live in libslatec.  Declare the dependency so the linker orders libslatec
# after libddeabm and resolves xermsg_ instead of dropping it.
target_link_libraries(ddeabm PUBLIC slatec)
