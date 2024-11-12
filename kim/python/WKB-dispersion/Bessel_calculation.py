from math import factorial
from scipy.special import gamma as Gamma
import numpy as np

def calc_needed_bessel_of_mphi(mphi, eval_b):
        
    BesselProd0 = 0 + 0j
    BesselProd1 = 0 + 0j
    BesselProd0_store = 1.0 + 0j
    BesselProd1_store = 1.0 + 0j
    k = 0
    while (np.abs(BesselProd0 - BesselProd0_store) > 1e-10 or np.abs(BesselProd1 - BesselProd1_store) > 1e-10):
        BesselProd0_store = 1.0 * BesselProd0
        BesselProd1_store = 1.0 * BesselProd1
        BesselProd0 += (0.5 * eval_b)**(2*k) / (factorial(k) * Gamma(k+mphi+1) + 0j) * np.exp(-eval_b)
        BesselProd1 += (0.5 * eval_b)**(2*k) / (factorial(k) * Gamma(k+mphi) + 0j) * np.exp(-eval_b)
        k +=1
        if k > 50:
            break
    BesselProd0 *= (0.5 * eval_b)**mphi
    BesselProd1 *= (0.5 * eval_b)**(mphi-1)
        
    return BesselProd0, BesselProd1 
    
#def calc_needed_bessel_of_mphi(self, mphi, eval_b):
    #if np.abs(np.real(eval_b))>self.bessel_large_arg_limit:
        ## asymptotic form for Bessel function:

        #BesselProd0 = 1.0 / (np.sqrt(2*np.pi*eval_b + 0j) * (1+ mphi**2 / eval_b**2)**(1/4)) * np.exp(- mphi * np.arcsinh(mphi / eval_b) + eval_b \
            #* np.sqrt(1 + mphi**2 / eval_b**2) - eval_b) + 0j
        #BesselProd1 = 1.0 / (np.sqrt(2*np.pi*eval_b + 0j) * (1+ (mphi-1)**2 / eval_b**2)**(1/4)) * np.exp(- (mphi-1) * np.arcsinh((mphi-1) / eval_b) + eval_b \
            #* np.sqrt(1 + (mphi-1)**2 / eval_b**2) - eval_b) + 0j
    #else:
        #BesselProd0 = Bessel(mphi, eval_b)* np.exp(-eval_b)
        #BesselProd1 = Bessel(mphi-1, eval_b)* np.exp(-eval_b)
    #return BesselProd0, BesselProd1 