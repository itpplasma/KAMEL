import sys

sys.path.append("../kilca_interface/")
from InpOut import InpOut


class field_divB0(InpOut):
    ###########################################################################
    # class field_divB0(InpOut)
    ###########################################################################
    # description of class:
    # --------------------------------------------------------------------------
    # This class is an interface for the field_divB0.f90 routine. It has the
    # content of the file field_divB0.inp (input file). Upon write() a new file
    # is created that can be used by the fortran routine. The method write()
    # is inherited from the InpOut class.
    ###########################################################################
    # properties:
    # --------------------------------------------------------------------------
    # *) ind, BLUEPRINT
    # *) ipert, iequil, ampl, ntor, cutoff, icftype
    # *) gfile, pfile, convexfile, fluxdatapath
    # *) winsize_Rfilt, winsize_Zfilt
    ###########################################################################
    # methods:
    # --------------------------------------------------------------------------
    # *) __init__(self, file_equi, file_coil, convex_path, fluxdp)
    # *) data(self) ... returns the data for the field_divB0.inp file
    ###########################################################################
    # More detailed comments on switches:
    # --------------------------------------------------------------------------
    # Axisymmetric equilibrium is called in all cases. If some routines are not
    # called, their data is not required.
    #
    # Perturbation switch "ipert":
    # 0 - axisymmetric equilibrium is called alone (can be run withot coild data
    #     and resonance plasma response data)
    # 1 - cylindrical routine is called in addition to equilibrium (plasma response
    #     data is not required
    # 2 - both, cylindrical and flux coordinate routine modelling plasma response
    #     are called, all data is required. Derivatives are not computed
    # 3 - the same as 2 with computation of derivatives (7 times slower)
    #
    # Equilibrium switch "iequil":
    # 0 - equilibrium field is not added to the output (needed for computing the
    #     perturbation field alone). The equilibrium routine is called - its
    #     data is needed for plasma response field. For future: this call can be
    #     disabled for stellarator modelling.
    # 1 - normal mode, perturbation is added to the output field
    ###########################################################################

    ind = list(range(0, 11 + 1))  # indices of parameters in blueprint files
    BLUEPRINT = "field_divB0.inp"
    # name of blueprint file
    sep = "  "

    ipert = 0  # 0=eq only, 1=vac, 2,3=vac+plas
    iequil = 1  # 0=pert. alone, 1=with equil.
    ampl = 1.00  # amplitude of perturbation, a.u.
    ntor = 72  # number of toroidal harmonics
    cutoff = 0.99  # inner cutoff in psi/psi_a units
    icftype = 4  # type of coil file
    gfile = ""  # equilibrium file (../3.0s/g33133.3000_ed4)
    pfile = ""  # field.dat - coil file (/itp/MooseFS/kasilov/MESH3D_31135odd/field.dat)
    convexfile = ""  # convex file for stretchcoords (./DATA/ASDEX/convexwall.dat)
    fluxdatapath = (
        ""  # directory with data in flux coord. (/proj/plasma/RMP/DATA2017/33133/4.3s/PREPROCESS)
    )
    winsize_Rfilt = 0  # window size for filtering of psi array over R
    winsize_Zfilt = 0  # window size for filtering of psi array over Z

    def __init__(self, file_equi, file_coil, convex_path, fluxdp):
        ###############################################################
        # function obj = field_divB0(gf, pf, convf, fluxdp)
        ###############################################################
        # description:
        # --------------------------------------------------------------
        # constructor of the field_divB0 class.
        ###############################################################
        # input:
        # --------------------------------------------------------------
        # gf    ... path to g-file (EFIT)
        # pf    ... path to p-file (field.dat)
        # convf ... path to convexfile
        # fluxdp... path to fluxdata directory
        ###############################################################

        self.gfile = file_equi
        self.pfile = file_coil
        self.convexfile = convex_path
        self.fluxdatapath = fluxdp

    def data(self):
        return [
            self.ipert,
            self.iequil,
            self.ampl,
            self.ntor,
            self.cutoff,
            self.icftype,
            self.gfile,
            self.pfile,
            self.convexfile,
            self.fluxdatapath,
            self.winsize_Rfilt,
            self.winsize_Zfilt,
        ]
