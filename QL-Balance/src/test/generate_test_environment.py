import numpy as np
from balance_interface import QL_Balance_interface
from utility import create_parabolic_profiles_from_res_surf

mpol = 6 # poloidal mode number
ntor = 2 # toroidal mode number
run_path = './test_run/'
profile_path = run_path + 'profiles/'

# all values at the resonant surface:
q0 = mpol/ntor # safety factor
n0 = 2e13 # cm^-3
Te0 = 1000 # eV
Ti0 = 1000 # eV
Vz0 = 1e6 # cm/s
Er0 = 0.2 # statV/cm
Vth0 = 1e5 # cm/s
rmin = 3.0
rmax = 80.0
num = 300 # number of radial grid points
a = 67.0 # plasma radius

Btor = -17000 # toroidal magnetic field on axis in Gauss

def main():
    bi = QL_Balance_interface(run_path = run_path, shot = 0, time = 0, name = 'test')
    bi.read_config_nml()
    bi.set_type_of_run(run_type = 'TimeEvolution')
    bi.set_modes(m_mode = mpol, n_mode = ntor)

    create_parabolic_profiles_from_res_surf(profile_path, q0, n0, Te0, Ti0, Vz0, Er0, Vth0, 
                              mpol, ntor, rmin, rmax, num, a)

    bi.prepare_balance(Btor=Btor, a_minor=a)
    bi.set_config_nml()
    bi.conf.conf['balancenml']['ramp_up_mode'] = 3 # instant max RMP coil current
    bi.conf.conf['balancenml']['t_max_ramp_up'] = 1.2
    bi.write_config_nml(path=run_path + 'balance_conf.nml')

if __name__ == "__main__":
    main()
