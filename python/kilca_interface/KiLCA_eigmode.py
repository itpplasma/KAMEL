from InpOut import InpOut
import numpy as np

class KiLCA_eigmode(InpOut):

    BLUEPRINT = 'eigmode.in'
    ind = ([2,5,22] + list(range(7,12+1)) + list(range(15,18+1)) + list(range(24, 25+1)) + list(range(28, 30+1)) + list(range(33,39+1))).sort()

    ouput = 'roots.dat'         # file name with determinant values
    flag_fscan = 1              # choose if to use frequency scan or roots search: 0 - if roots search is used, 1 - if det is evaluated on a frequency grid

    # frequency grid settings:
    fgrid_redim = 1             # dimension of a grid for real part
    fgrid_remin = 0             # minimum real value of frequency
    fgrid_remax = 1e6           # maximum real value of frequency
    fgrid_imdim = 100           # dimension of a grid for imag part
    fgrid_immin = 0             # minimum imag value of frequency
    fgrid_immax = 2.1e5         # maximum imag value of frequency

    #stopping criteria for roots search:
    flag_stopcrit = 1           # stopping criterion: 0 -det residual, 1 - root convergence
    det_abs = 1e-14             # absolute value of the determinant
    rootseq_abserr = 1e-10      # absolute error for roots sequence
    rootseq_relerr = 1e-10      # relative error for roots sequence

    # for omega derivative:
    df = 1e-3                   # delta freq used to compute derivative over omega numerically

    # test roots:
    flag_testroots = 0          # 1 - test roots by contour integration, 0 - skip test
    flag_deb = 0                # flag for debugging

    # starting point indices for roots search
    rsearch_nstart = 4          # number of starting points
    rsearch_istart = 0          # starting index (from 0)
    rsearch_iend = 3            # ending index (from 0)
    rsearch_points = np.array((np.zeros(7), np.logspace(2,8,7))).T.tolist()

    def __init__(self):
        pass

    def data(self):
        return [self.ouput, self.flag_fscan, self.fgrid_redim, self.fgrid_remin, self.fgrid_remax, self.fgrid_imdim, self.fgrid_immin, self.fgrid_immax, self.flag_stopcrit, self.det_abs, self.rootseq_abserr, self.rootseq_relerr, self.df, self.flag_testroots, self.flag_deb, self.rsearch_nstart, self.rsearch_istart, self.rsearch_iend, self.rsearch_points]
