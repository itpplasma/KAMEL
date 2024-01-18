# Analytical version of the approximate bifurcation criterion (A2BC)
# The approximate bifurcation criterion is given by D^{ql}_{e22}/D^{a}>1.
# The function a2bc returns the ratio of the diffusion coefficients for a given set of input data.
# author: Markus Markl
# date: 02.01.2023

import numpy as np
from scipy.interpolate import CubicSpline

sol = 29979245800 # speed of light
const_e = 4.8032e-10 # elementary charge
const_kB = 1.38007e-16 # Boltzmann constant
const_eVK = 1.604e4 # conversion factor eV to Kelvin
const_eV_to_erg = 1.6022e-12
e_mass = 9.1094e-28#511

debug = False

def analytical_local_criterion(m, n, r_Da, Da, R0, B0, r, q, Te, Ti, ne, Er, Impar, Zi=1, mi=2):
    '''Analytical local bifurcation criterion.'''
    # all quantities need to be in cgs units, except for the temperatures
    # which need to be in eV
    # m ... poloidal mode number (int)
    # n ... toroidal mode number (int)
    # r_Da ... radial variable of the anomalous diffusion coefficient (array)
    # Da ... anomalous diffusion coefficient (array or double)
    # R0 ... major radius of the device (double)
    # r ... radial coordinate of profiles (array)
    # q ... safety factor profile (array or double)
    # Te ... electron temperature profile (array or double)
    # ne ... electron density profile (array or double)
    # Er ... radial electric field profile (array or double)
    # Impar ... parallel plasma response current from e.g. GPEC or KiLCA (double)
    # Zi ... charge number ions, default=1, (int)
    # mi ... mass number ions, default=2, (int)

    if type(r) == np.ndarray:
        rm = get_rm(r,m,n,q)
    else:
        rm = r

    B0 = np.abs(B0)

    # Check if input is array or scalar
    if type(Te) == np.ndarray:
        Te_res = np.interp(rm, r, Te)
    else:
        Te_res = Te
    if type(Ti) == np.ndarray:
        Ti_res = np.interp(rm, r, Ti)
    else:
        Ti_res = Ti
    if type(ne) == np.ndarray:
        ne_res = np.interp(rm, r, ne)
    else:
        ne_res = ne
    if type(Er) == np.ndarray:
        Er_res = np.interp(rm, r, Er)
    else:
        Er_res = Er
    if type(q) == np.ndarray:
        q_res  = np.interp(rm, r, q)
    else:
        q_res = q
    if type(Da) == np.ndarray:
        Da_res = np.interp(rm, r_Da, Da)
    else:
        Da_res = Da

    #rDe = 7.43e2 * Te_res**(1/2) * ne_res**(-1/2) # from NRL formulary 2019, is less accurate
    rDe = np.sqrt(Te_res * const_eV_to_erg/(4 * np.pi * ne_res * const_e**2)) # Debye length in cm
    kz = n / R0
    kperp = -(m / rm - (n * rm/ (R0**2 * q_res))) # Bz/B0 \approx -1
    omega_E = - kperp * sol * Er_res / B0
    
    # magnetic shear
    s_res = rm / q_res * np.interp(rm, r, np.gradient(q,r))
    
    #Coulomb logarithms from KiLCA (NRL formulary 2009)
    Lee = 23.5 - np.log(np.sqrt(ne_res) / Te_res**1.25) - np.sqrt(1e-5 + (np.log(Te_res) - 2.0)**2 / 16)
    Lei = 24 - np.log(np.sqrt(ne_res) *Ti_res**(-1.0))
    
    # collision frequencies
    nuee = 5.8e-6 * ne_res * Lee/Te_res**(3/2)
    nuei = 7.7e-6 * ne_res * Lei * Zi**2 / Te_res**(3/2)
    nue = nuee + nuei + 10

    # diamagnetic velocity of electrons
    VeD = sol / (const_e * ne_res * B0) * np.interp(rm, r, np.gradient(ne * Te * const_eV_to_erg, r))
    # ExB velocity
    VE = sol * Er_res / B0
    # temperature part of the diamagnetic velocity
    VDTe = sol / (const_e * B0) * np.interp(rm, r, np.gradient(Te* const_eV_to_erg, r))
    # resonance combination of the velocities
    Veres = VeD + VE - 1/2 * VDTe
    #print('other = ' + str(0.5 * VDTe - VeD))
    #print('VE    = ' + str(VE))
    #print('Veres = ' + str(Veres))

    # parts of the bifurcation criterion
    crit1 = (rDe / rm)**4
    crit2 = s_res**2 * sol**2 / (8 * np.pi**2 * Veres**2)
    crit3 = np.abs((2 * np.pi * omega_E - 6 * 1j* nue) / omega_E)
    crit4 = (279 * nue**3 + 47 * omega_E**2 * nue) / (9 * nue**4 + 10 * omega_E**2 * nue**2 + omega_E**4)
    crit5 = np.abs(Impar * sol)**2 * kz**2 / (Da_res * B0**2)
    crit = crit1*crit2*crit3*crit4*crit5 # factor 2 is a spectral weight since the Fourier decomposition of the fields in KiLCA are only over positive integers, but GPEC uses pos and negative integers
    

    
    if debug:
        print('a2bc')
        print('rm      = ' + str(rm))
        print('Te res  = ' + str(Te_res))
        print('Ti res  = ' + str(Ti_res))
        print('ne res  = ' + str(ne_res))
        print('Er res  = ' + str(Er_res))
        print('q res   = ' + str(q_res))
        print('Da res  = ' + str(Da_res))
        print('rDe     = ' + str(rDe))
        print('kz      = ' + str(kz))
        print('kperp   = ' + str(kperp))
        print('omega E = ' + str(omega_E))
        print('s_res   = ' + str(s_res))
        print('Lee     = ' + str(Lee))
        print('Lei     = ' + str(Lei))
        print('nue     = ' + str(nue))
#        print('omega E / nue = ' + str(omega_E / nue))
        print('VeD     = ' + str(VeD))
        print('VE      = ' + str(VE))
        print('VDTe    = ' + str(VDTe))
        print('Veres   = ' + str(Veres))
        print('crit1   = ' + str(crit1))
        print('crit2   = ' + str(crit2))
        print('crit3   = ' + str(crit3))
        print('crit4   = ' + str(crit4))
        print('crit5   = ' + str(crit5))
        print('crit    = ' + str(crit))

    return crit


def Br_const_psi(m, n, R0, B0, r, q, Te, Ti, ne, Er, Impar, Zi=1, mi=2):
    '''Calculate Br in constant psi approimation.'''

    rm = get_rm(r,m,n,q)

    B0 = np.abs(B0)

    Te_res = np.interp(rm, r, Te)
    Ti_res = np.interp(rm, r, Ti)
    ne_res = np.interp(rm, r, ne)
    Er_res = np.interp(rm, r, Er)
    q_res  = np.interp(rm, r, q)

    #rDe = 7.43e2 * Te_res**(1/2) * ne_res**(-1/2) # from NRL formulary 2019, is less accurate
    rDe = np.sqrt(Te_res * const_eV_to_erg/(4 * np.pi * ne_res * const_e**2)) # Debye length in cm
    kz = n / R0
    kperp = -(m / rm - (n * rm/ (R0**2 * q_res))) # Bz/B0 \approx -1
    omega_E = - kperp * sol * Er_res / B0
    
    # magnetic shear
    s_res = rm / q_res * np.interp(rm, r, np.gradient(q,r))
    
    #Coulomb logarithms from KiLCA (NRL formulary 2009)
    Lee = 23.5 - np.log(np.sqrt(ne_res) / Te_res**1.25) - (np.log(Te_res) - 2.0)/4.0
    Lei = 30.0 - np.log(np.sqrt(ne_res) / Ti_res**1.5 * Zi**2 / mi)
    
    # collision frequencies
    nuee = 5.8e-6 * ne_res * Lee/Te_res**(3/2)
    nuei = 7.7e-6 * ne_res * Lei * Zi**2 / Te_res**(3/2)
    nue = nuee + nuei + 10

    # diamagnetic velocity of electrons
    VeD = sol / (const_e * ne_res * B0) * np.interp(rm, r, np.gradient(ne * Te * const_eV_to_erg, r))
    # ExB velocity
    VE = sol * Er_res / B0
    # temperature part of the diamagnetic velocity
    VDTe = sol / (const_e * B0) * np.interp(rm, r, np.gradient(Te* const_eV_to_erg, r))
    # resonance combination of the velocities
    Veres = VeD + VE - 1/2 * VDTe

    #vTe = 4.19e7 * Te_res**(1/2)
    vTe = np.sqrt(Te_res  * const_eV_to_erg/ e_mass)

    bcp1 = - (s_res * kz * sol * rDe**2) / (np.pi * vTe * rm**2)
    bcp2 = np.sqrt(2 * np.pi - 6 * 1j * nue/omega_E)
    bcp3 = Impar * sol / Veres
    bcp = bcp1 * bcp2 * bcp3


    #debug = False

    if debug: print('Br_const_psi')
    if debug:
        print('rm      = ' + str(rm))
        print('Te res  = ' + str(Te_res))
        print('Ti res  = ' + str(Ti_res))
        print('ne res  = ' + str(ne_res))
        print('Er res  = ' + str(Er_res))
        print('q res   = ' + str(q_res))
        print('rDe     = ' + str(rDe))
        print('kz      = ' + str(kz))
        print('kperp   = ' + str(kperp))
        print('omega E = ' + str(omega_E))
        print('s_res   = ' + str(s_res))
        print('Lee     = ' + str(Lee))
        print('Lei     = ' + str(Lei))
        print('nue     = ' + str(nue))
#        print('omega E / nue = ' + str(omega_E / nue))
        print('VeD     = ' + str(VeD))
#        print('    VeD at res surf.     = ' + str(np.interp(rm, r, VeD)))
        print('VE      = ' + str(VE))
#        print('    VE at res surf.     = ' + str(np.interp(rm, r, VE)))
        print('VDTe    = ' + str(VDTe))
#        print('    VDTe at res surf.     = ' + str(np.interp(rm, r, VDTe)))
        print('Veres   = ' + str(Veres))
#        print('    Veres at res surf.     = ' + str(np.interp(rm, r, Veres)))
        print('vTe     = ' + str(vTe))
        print(' - - - ')

    return bcp


def crit_from_Br(m, n,r_Da, Da, R0, B0, r, q, Te, Ti, ne, Er, Impar, Br = 0, Zi=1, mi=2):
    '''Determine criterion from a given Br.'''

    if Br==0:
        bcp = Br_const_psi(m, n, R0, B0, r, q, Te, Ti, ne, Er, Impar, Zi=Zi, mi=mi)
    else:
        bcp = Br

    rm = get_rm(r,m,n,q)

    #Te_res = np.interp(rm, r, Te)
    #Ti_res = np.interp(rm, r, Ti)
    #ne_res = np.interp(rm, r, ne)
    #Er_res = np.interp(rm, r, Er)
    #q_res  = np.interp(rm, r, q)
    Da_res = np.interp(rm, r_Da, Da)

    kperp = -(m / rm - (n * rm/ (R0**2 * q))) # Bz/B0 \approx -1
    omega_E = - kperp * sol * Er / B0
    
    #Coulomb logarithms from KiLCA (NRL formulary 2009)
    Lee = 23.5 - np.log(np.sqrt(ne) / Te**1.25) - (np.log(Te) - 2.0)/4.0
    Lei = 30.0 - np.log(np.sqrt(ne) / Ti**1.5 * Zi**2 / mi)
    
    # collision frequencies
    nuee = 5.8e-6 * ne * Lee/Te**(3/2)
    nuei = 7.7e-6 * ne * Lei * Zi**2 / Te**(3/2)
    nue = nuee + nuei + 10

    #vTe = 4.19e7 * Te_res**(1/2)#np.sqrt(Te_res  * const_kB * const_eVK/ e_mass)

    vTe = np.sqrt(Te  * const_eV_to_erg/ e_mass)

    crit1 = vTe**2 / 8
    crit2 = (47 * omega_E**2 * nue + 279*nue**3)/(omega_E**4 + 10 * omega_E**2 * nue**2 + 9 * nue**4)
    crit3 = np.abs(bcp)**2 / B0**2

    crit = crit1*crit2*crit3 / Da_res
    crit = np.interp(rm, r, crit)

    #debug = False

    if debug: print('crit_from_Br')
    if debug:
        print('values at res surf:')
        print('rm      = ' + str(rm))
        print('kperp   = ' + str(np.interp(rm, r, kperp)))
        print('omega E = ' + str(np.interp(rm, r, omega_E)))
        print('Lee     = ' + str(np.interp(rm, r, Lee)))
        print('Lei     = ' + str(np.interp(rm, r, Lei)))
        print('nue     = ' + str(np.interp(rm, r, nue)))
        print('vTe     = ' + str(np.interp(rm, r, vTe)))
        print(' - - - ')

    return crit

def get_rm(r,m,n,q):
    '''Get the radius of a rational surface between 45 < r < 65 cm.'''

    ind2 = np.where(r< 65)

    q_new = q[ind2]
    r_out_new = r[ind2]

    ind1 = np.where(r_out_new> 45)
    q_new = q_new[ind1]
    r_out_new = r_out_new[ind1]

    f = CubicSpline(np.abs(q_new), r_out_new)

    return f(m/n)