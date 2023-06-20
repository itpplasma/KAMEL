from InpOut import InpOut
from KiLCA_zone_vacuum import KiLCA_zone_vacuum
from KiLCA_zone_flre import KiLCA_zone_flre

class KiLCA_zone(InpOut):
    """
    Description of the class:
        Class containing the information of a zone in KiLCA. Corresponds to a zoneX.in
        file for KiLCA.
    ############################
    Variables:
        r1, typeBC1, model, modelvers, typeBC2, r2
        vacuum, imhd, flre
        ind, BLUEPRINT, READY
    Methods:
        __init__(num, r1, b1, m, r2, b2)
        plain()
        indices
    """

    BLUEPRINT = ''
    ind = []
    sep = '#'

    number = 0 # number of zones
    
    r1 = []      # r1 - minimum radius of the zone (plasma radius)
    typeBC1 = [] # type of BC at r1 (center, infinity, interface, antenna, idealwall)
    model = []   # type of the plasma model (vacuum, medium, imhd, rmhd, flre)
    modelvers = 0# code version for the model: MHD model (0- incompressible and flowless, 1 - compressible with flows)
    typeBC2 = [] # type of BC at r2 (center, infinity, interface, antenna, idealwall)
    r2 = []      # r2 - maximum radius of the zone (first wall boundary)

    vacuum = []  # contains vacuum information about the zone
    imhd   = []  # contains imhd information about the zone
    flre   = []  # contains flre information about the zone

    typeBC1_num = 0 # numeric type of BC at r1 (center=0, infinity=1, interface=2, antenna=3, idealwall=4)
    model_num = 0   # numeric type of the plasma model (vacuum=0, medium=1, imhd=2, rmhd=3, flre=4)
    typeBC2_num = 0 # numeric type of BC at r2 (center=0, infinity=1, interface=2, antenna=3, idealwall=4)



    def __init__(self, num, r1, b1, m, r2, b2):
        """
        Constructor of the KiLCA_zone class.
        Input:
            num... number of zones
            r1 ... inner radius
            b1 ... inner boundary type
            m  ... model between borders (vacuum, flre, mhd)
            r2 ... outer radius
            b2 ... outer boundary type
        """
        self.r1 = r1
        self.typeBC1 = b1
        self.model = m
        self.r2 = r2
        self.typeBC2 = b2

        if m =='imhd' or m=='rmhd':
            raise ValueError('imhd or rmhd support not yet available')

        if m =='vacuum':
            self.vacuum = KiLCA_zone_vacuum()
            o = self.vacuum
        elif m=='flre':
            self.flre = KiLCA_zone_flre()
            o = self.flre
        
        self.number = num
        self.ind = o.ind
        self.BLUEPRINT = o.BLUEPRINT

        self.READY = True

        

    def get_typeBC1_num(self):
        return self.BC_num(self.typeBC1)
    def get_typeBC2_num(self):
        return self.BC_num(self.typeBC2)
    def get_model_num(self):
        return self.MD_num(self.model)

    def BC_num(self, prop):
        """
        Boundary condition to number.
        """
        if prop == 'center':
            return 0
        elif prop == 'infinity':
            return 1
        elif prop == 'interface':
            return 2
        elif prop == 'antenna':
            return 3
        elif prop == 'idealwall':
            return 4

    def MD_num(self, prop):
        """
        Model numeric type.
        """
        if prop == 'vacuum':
            return 0
        elif prop == 'medium':
            return 1
        elif prop == 'imhd':
            return 2
        elif prop == 'rmhd':
            return 3
        elif prop == 'flre':
            return 4