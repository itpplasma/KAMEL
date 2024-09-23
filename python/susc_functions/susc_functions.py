import numpy as np
import susc_funcs
import matplotlib.pyplot as plt
from plasmapy.dispersion import plasma_dispersion_func as plasma_disp


x1 = np.linspace(-50,50,10000)
x2 = np.logspace(-2,1,5)

if False:
    plt.figure()
    for x_2 in x2:
        I00 = np.zeros(len(x1))
        for i, x in enumerate(x1):
            I = susc_funcs.getifunc(np.abs(x),x_2)
            I00[i] = np.real(I[0,0])
        plt.plot(np.abs(x1),np.real(I00)/np.max(np.real(I00)), label='x2 = %.2f' % x_2)
    plt.legend(bbox_to_anchor=(1.05,1))
    plt.xlabel('x1')
    plt.ylabel('Re(I00)/max(Re(I00))')
    plt.tight_layout()
    plt.show()

Imn_arr = susc_funcs.calc_imn_array(x1,x2)
print(Imn_arr)

x_2 = 20

def compare_I20_appro():
    plt.figure()
    I00 = np.zeros(len(x1), dtype=complex)
    I00_arr = np.zeros(len(x1), dtype=complex)
    I20_approx = np.zeros(len(x1), dtype=complex)
    for i, x in enumerate(x1):
        I = susc_funcs.getifunc(x,x_2)
        I00[i] = np.real(I[2,0])
        I = susc_funcs.calc_imn_array(x,x_2)
        I00_arr[i] = np.real(I[2,0])
        I20_approx[i] = - 1j * x_2 /x**2 * (x_2/(np.sqrt(2) * x) * np.sign(x) * plasma_disp(x_2/(np.sqrt(2) * x)) + 1.0)

    plt.plot(x1, np.real(I00), label='I20 (Balance)')
    plt.plot(x1, np.real(I00_arr), label='I20 arr (KiLCA)')
    plt.plot(x1, np.real(I20_approx), label='I20 aapprox')

    plt.legend(bbox_to_anchor=(1.05,1))
    plt.show()

def compare_I00_appro():
    plt.figure()
    I00 = np.zeros(len(x1), dtype=complex)
    I00_arr = np.zeros(len(x1), dtype=complex)
    I00_approx = np.zeros(len(x1), dtype=complex)
    for i, x in enumerate(x1):
        I = susc_funcs.getifunc(x,x_2)
        I00[i] = np.real(I[0,0])
        I = susc_funcs.calc_imn_array(x,x_2)
        I00_arr[i] = np.real(I[0,0])
        I00_approx[i] = -1j * np.sign(x) /(np.sqrt(2) * x) * plasma_disp(x_2/(np.sqrt(2) * x))

    plt.plot(x1, np.real(I00), label='I00 (Balance)')
    plt.plot(x1, np.real(I00_arr), label='I00 arr (KiLCA)')
    plt.plot(x1, np.real(I00_approx), label='I00 aapprox')

    plt.legend(bbox_to_anchor=(1.05,1))
    plt.tight_layout()
    plt.show()

compare_I00_appro()