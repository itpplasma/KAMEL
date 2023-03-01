from InpOut import InpOut

class KiLCA_background(InpOut):
    """
    This class handles the background configuration.
    """

    BLUEPRINT = 'background.in'
    ind = list(range(1,3+1)) + list(range(6,14+1)) + [17]
    sep = '#'

    Rtor = 165.0        # big torus radius (cm) of the machine
    rpl = 63.0          # plasma radius (cm)
    Btor = 2e4          # toroidalmagnetic field (G) at the center
    prof_path = './profiles/' # path to input background profiles
    flag_recalc = 1     # 1 - if background must be recalculated (7 input profiles are needed), 0 - otherwise (11 input profiles are needed), -1 -from interface
    flag_back = 'f'     # flag for background ('f'-full, 'w'-wkb, 'h'-hom)
    splinedeg = 9       # splines degree: >= NC + 2N +1, where N - order of FLR expansion, NC - spl degree for C matrices, must be odd!!
    vgalsys = -1e9      # V_gal_sys is a velocity (cm/s) of a moving frame
    vscale = 1          # V_scale: scale factor for the Vz velocity profile: Vz = V_scale * Vz - V_gal_sys
    mi = 2              # m_i: ions mass in units of proton mass
    ce = 1              # collisions coefficient for electrons
    ci = 1              # collisions coefficient for ions
    flag_deb = 0        # flag for debugging mode (additional checks are performed in the code)

    def __init__(self, R, r):
        """
        Constructor takes big torus radius and plasma radius.
        """

        self.Rtor = R
        self.rpl = r

    def data(self):
        return [self.Rtor, self.rpl, self.Btor, self.prof_path, self.flag_recalc, self.flag_back, self.splinedeg, self.vgalsys, self.vscale, self.mi, self.ce, self.ci, self.flag_deb]