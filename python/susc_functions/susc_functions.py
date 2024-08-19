import numpy as np
import susc_funcs
import matplotlib.pyplot as plt


x1 = np.linspace(-10,10,10000)
x2 = np.logspace(-2,1,5)

plt.figure()
for x_2 in x2:
    I00 = np.zeros(len(x1))
    for i, x in enumerate(x1):
        I = susc_funcs.getifunc(x,x_2)
        I00[i] = np.real(I[0,0])
    plt.plot(x1,np.real(I00)/np.max(np.real(I00)), label='x2 = %.2f' % x_2)
plt.legend(bbox_to_anchor=(1.05,1))
plt.xlabel('x1')
plt.ylabel('Re(I00)/max(Re(I00))')
plt.tight_layout()
plt.show()