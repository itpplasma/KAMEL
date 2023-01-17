from InpOut import InpOut

class KiLCA_antenna(InpOut):
    """
    This class handles the antenna data. It is used to read and write the antenna.in text file used in KiLCA.
    """

    BLUEPRINT = 'antenna.in'
    ra          = 70         # radius (cm) of the antenna location (must match zones.in)
    width       = 0.0     # current density layer width, if ==0 then use delta function
    I0          = 4.5e12 # current in the coils (statamp). 4.5e12 same as 1.5kA
    flab = [1e3,0] # complex frequency (1/c) in the laboratory frame
    nmod        = 1    # number of antenna modes to be used from modes.in file
    flag_deb    = 0 # flag for debugging
    flag_eig    = 0 # flag to solve eigenmode problem

    def __init__(self, r, n):
        """
        Constructor takes radius r of antenna and number of modes n.
        """

        self.ra = r
        self.nmod = n
        self.READY = True

    def export_to_hdf5(self, fname, loc):
        pass

    def export_to_nml(self, fname, loc):
        pass
