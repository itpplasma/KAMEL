from InpOut import InpOut
from KiLCA_zone_vacuum import KiLCA_zone_vacuum
from KiLCA_zone_flre import KiLCA_zone_flre
from read_in import read_in
from change_opts import change_opts
from save_file import save_file

class KiLCA_zone(InpOut):
    """
    Description of the class:
        Class containing the information of a zone in KiLCA. Corresponds to a zone_*.in
        file for KiLCA. Inherits InpOut class to read and write files.
    ############################
    Variables:
        boundary_cond, model_type
        ind, BLUEPRINT, READY
    Methods:
        __init__(num, r1, b1, m, r2, b2)
            Constructor of the class. Writes arguments into data dict and creates appropriate
            KiLCA zone instance of vacuum or flre class.
        get_typeBC1_num()
            Returns index of boundary condition 1.
        get_typeBC2_num()
            Returns index of boundary condition 2.
        get_model_num()
            Returns index of model.
        return_data()
            Method to overwrite the inherited method from InpOut since data is a combination
            of zone description (r1, typeBC1,...) as well as zone class data.
        return_description_of_zone()
            Invoked in return data. Returns the description of the zone (r1, typeBC1,...)
    """

    BLUEPRINT = ''
    ind = []
    sep = '#'

    boundary_cond = ['center', 'infinity', 'interface', 'antenna', 'idealwall']
    model_type = ['vacuum', 'medium', 'imhd', 'rmhd', 'flre']
    o = None

    def __init__(self, num: int, r1: float, b1: str, m: str, r2: float, b2: str):
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

        self.data = {
            'number': 0, # number of zones
            'r1': [], # r1 - minimum radius of the zone (plasma radius)
            'typeBC1': [], # type of BC at r1 (center, infinity, interface, antenna, idealwall)
            'model': [],   # type of the plasma model (vacuum, medium, imhd, rmhd, flre)
            'modelvers': 0, # code version for the model: MHD model (0- incompressible and flowless, 1 - compressible with flows)
            'typeBC2': [], # type of BC at r2 (center, infinity, interface, antenna, idealwall)
            'r2': [],      # r2 - maximum radius of the zone (first wall boundary)
            #'vacuum': [], # contains vacuum information about the zone
            #'imhd': [],  # contains imhd information about the zone
            #'flre': [],  # contains flre information about the zone

            'typeBC1_num': 0, # numeric type of BC at r1 (center=0, infinity=1, interface=2, antenna=3, idealwall=4)
            'model_num': 0,   # numeric type of the plasma model (vacuum=0, medium=1, imhd=2, rmhd=3, flre=4)
            'typeBC2_num': 0 # numeric type of BC at r2 (center=0, infinity=1, interface=2, antenna=3, idealwall=4)
        }
        self.data['r1'] = r1
        self.data['typeBC1'] = b1
        self.data['model'] = m
        self.data['r2'] = r2
        self.data['typeBC2'] = b2

        if m =='imhd' or m=='rmhd':
            raise ValueError('imhd or rmhd support not yet available')
        elif m =='vacuum':
            self.data['vacuum'] = KiLCA_zone_vacuum()
            self.o = self.data['vacuum']
        elif m=='flre':
            self.data['flre'] = KiLCA_zone_flre()
            self.o = self.data['flre']
        
        self.data['number'] = num
        self.ind = self.o.ind
        self.BLUEPRINT = self.o.BLUEPRINT

        self.READY = True


    def get_typeBC1_num(self):
        return self.boundary_cond.index(self.data['typeBC1'])

    def get_typeBC2_num(self):
        return self.boundary_cond.index(self.data['typeBC2'])

    def get_model_num(self):
        return self.model_type.index(self.data['model'])
        
    def write(self, path_from, path_to):
        """
        Description:
            writes propeties of the class into input file.
            Overrides superclass method because output file name is different from file name.
        Input:
            path_from ... path where the blueprint is from
            path_to   ... path where the input file will be written
        """

        if self.READY == False:
            raise ValueError('Class is not ready to run')

        # read blueprint

        raw = read_in(path_from + self.BLUEPRINT)
        raw = change_opts(raw, self.ind, self.return_data(), self.sep)
        save_file(raw, path_to + '/' + (self.BLUEPRINT).replace('zone', 'zone_'+str(int(self.data['number'])+1)))

    def return_data(self):
        return self.data_description_of_zone() + self.o.return_data()

    def data_description_of_zone(self):
        return [self.data['r1'], self.data['typeBC1'], self.data['model'], self.data['modelvers'], self.data['typeBC2'], self.data['r2']]
