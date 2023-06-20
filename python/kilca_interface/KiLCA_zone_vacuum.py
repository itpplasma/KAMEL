from InpOut import InpOut

class KiLCA_zone_vacuum(InpOut):
    """
    Description:
        Class representing a vacuum zone in KiLCA.

    Variables:
        sigma
        rgrid_maxdim, relacc, absacc
        polydeg, sparse_relacc, sparse_absacc, maxgridstep
        flag_deb
        ind, BLUEPRINT, READY

    Methods:
        KiLCA_zone_vacuum()
        plain()
    """

    ind = (list(range(1,6+1)) + [9,23] + list(range(12,14+1)), list(range(17, 20+1))).sort()
    BLUEPRINT = 'zone_vacuum.in'
    SEP = '#'

    sigma           = [0, 0] # conductivity
    rigid_maxdim    = 1000   # max dimension of the radial grid for the solution: default=1e5
    relacc          = 1e-8   # relative accuracy of the solution: default=1e-8
    absacc          = 1e-8   # absolute accuracy of the solution: default=1e-8

    # ME solution space out settings
    polydeg         = 3      # degree of the polynomial used to space out the solution (by checking the accuracy)
    sparse_relacc   = 1e-8   # relative accuracy of the sparse solution: default=1e-8
    sparse_abs_acc  = 1e-8   # absolute accuracy of the sparse solution: default=1e-8
    maxgridstep     = 0.1    # max grid step in the solution: default=0.1

    flag_deb        = 0      # flag for debugging mode


    def __init__(self):
        """
        Description:
            Constructor of the vacuum zone class. Empty.
        """
        pass