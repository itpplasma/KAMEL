from InpOut import InpOut

class KiLCA_output(InpOut):
    """
    Description:
        Class containing the information of the output in KiLCA.
        General run settings: 0 - to skip, 1 - to calculate only, 2 - to calculate and save to disk
    Variables:
        data... contains all relevant data for the output.in file
        ind, BLUEPRINT, sep
    Methods:
        KiLCA_output()
    Comments:
        If appropriate the quantity is computed for ions (i), electrons (e) and total (i+e).
        If appropriate the quantity is computed for several types of the current, j0 (0) - basic current density
        Bear in mind that some quantities are dependent on others.
    """

    ind = sorted(list(range(1,4+1)) + [7] + list(range(10,17+1)))
    BLUEPRINT = 'output.in'
    sep = '#'

    data = {
        'backdata': 2,     #background data: default=2
        'lindata': 2,     #linear data: default=2
        'varquant': 2,     #various quantities (see below): default=2
        'disp': 0,     #dispersion: default=0
        'flag_deb': 0,     #flag for debugging
        'curdenspet': 1,     #flag: current density perturbation
        'abspow': 1,     #flag: absorbed power (density and integrated over cylinder volume)
        'dispow': 1,     #flag: dissipated power (density and integrated over cylinder volume)
        'kinfluxr': 1,     #flag: r-component of kinetic flux out of the cylinder volume
        'poyntfluxr': 1,     #flag: r-component of Poynting flux out of the cylinder volume
        'totfluxr': 1,     #flag: r-component of total flux out of the cylinder volume
        'numdenspet': 1,     #flag: number density perturbation
        'lortorq': 1     #flag: Lorentz torque (density and integrated over cylinder volume)
    }

    def __init__(self):
        pass
