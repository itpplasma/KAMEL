from .InpOut import InpOut

class KiLCA_zone_flre(InpOut):
    """
    Description:
        Class representing a flre zone in KiLCA.

    Variables:
        sigma
        rgrid_maxdim, relacc, absacc
        polydeg, sparse_relacc, sparse_absacc, maxgridstep
        flag_deb
        ind, BLUEPRINT, READY
    
    Methods:
        KiLCA_zone_flre()
        plain(self)

    Comments:
        If the flag for debugging mode is 1 then additional checks are performed in the code
        if it is 2 then additional checks are made, if 0 - no checks.

        Controlling factor for ONS is rather important parameter, default values are good only for the N=1 order, for higher orders of FLRE it might be necessary to decrease it up to 3.0 or so and also increase accuracy of integration up to 1.0e-12. The quality of the solution can be checked by inspection of bp.dat file - the solution must be SMOOTH! The file poy_test_err.dat is also useful
    """

    ind = sorted(list(range(1,6+1)) + list(range(9, 17+1)) + list(range(20, 27+1)) + list(range(30, 33+1)) + [36] + [39,40])
    BLUEPRINT = 'zone_flre.in'

    data = {
        'order': 1, # order of the FLR expansion: must match to sources used for conductivity
        'max_cycharm': 1, # highest cyclotron harmonic normally should exceed the flre order
        'flag_corr': 1, # flag if correction term should be used in conductivity: use for flre order =1 only
        'splinedeg': 5, # splines degree for C matrices, the degree for K matrices is computed in the code: NK = NC + N + 1
        'cond_maxdim': 3001, # maximum dimension of the radial grid for conductivity matrices: default=3001
        'reswidth': 3.0, # resonant layer width
        'outer_err': 1e-6, # error parameter used for adaptive radial grid generation outside the resonant layer: default=1.0e-6
        'inner_err': 1e-6, # error parameter used for adaptive radial grid generation in the resonant layer: default=1.0e-6
        'flag_seqhom': 0, # flag if system of equations should be used in homogeneous limit: for flre order greater than 1
        'ode_maxdim': 1e5, # max dimension of the radial grid for the solution: default=1e5
        'ode_relacc': 1e-8, # relative accuracy of the solution by ODE solver: default=1e-8
        'ode_absacc': 1e-8, # absolute accuracy of the solution by ODE solver: default=1e-8
        'ode_nons': 50000, # max number of the orthonormalization steps (ONS) for the solver: default=50000,
         # needs to be written like this (integer)! the exponential form yields a float
        'ode_cont': 1e3, # controlling factor for ONS by QR: norm_max/norm_min > norm_fac: default =1e3
        'ode_outer_step': 1e-3, # output grid step outside the resonance region for the ME solutions: default=1.0e-3
        'ode_inner_step': 1e-5, # output grid step inside the resonance region for the ME solutions: default=1.0e-5
        'ode_reswidth': 1, # width of the resonance region: default = 1.0
        'polydeg': 5, # degree of the polynomial used to space out the solution (by checking the accuracy)
        'sparse_relacc': 1e-8, # relative accuracy of the sparse solution: default = 1e-8
        'sparse_absacc': 1e-8, # absolute accuracy of the sparse solution: default = 1e-8
        'maxgridstep': 0.1, # max grid step in the solution: default=0.1
        'flag_deb': 0, # flag for debugging mode (additional checks in the code are done)
        'i_col': 1, # collisions model for ions: 0 -(N) = cont, 1 - (n,E) = const, 2 - (N,P) = const, 3 - (N, P,E) = const, where N - number, E - energy, P - momentum
        'e_col': 1 # collisions model for electrons: 0 -(N) = cont, 1 - (n,E) = const, 2 - (N,P) = const, 3 - (N, P,E) = const, where N - number, E - energy, P - momentum
    }

    def __init__(self) -> None:
        """
        Description:
            Constructor of the flre zone class. Empty.
        """
        pass

