class AUG_config:

    R0 = 165.0
    r_eff_wall = 80.0  # ideal wall effective radius
    r_eff_plas = 67.0
    delta_r_antenna = 3.0  # distance of RMP antenna from plasma boundary in cm
    r_eff_antenna = r_eff_plas + delta_r_antenna  # antenna effective radius
    Btor = -17563.3704  # shot dependent, but gives order of magnitude
    I0_rmp = 5.1e12

    def __init__(self):
        pass
