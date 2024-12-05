# terrible style - NO GLOBAL VARIABLES!

import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

N = 2024
Dx = 0.05
Lx = N*Dx
Dxi = 2*np.pi/Lx
Lxi = N*Dxi
x = Dx*np.arange(-N/2, N/2)
xi = Dxi*np.arange(-N/2, N/2)


# Parameters
k = 5  # degrees of freedom
lambda_ = 5  # non-centrality parameter

Fa = np.divide(np.exp(np.divide(1j*lambda_*xi,1-2j*xi)), pow(1-2j*xi,k/2))

fn = fftshift(fft(ifftshift(Fa)))/Lx


plt.figure(1)
plt.plot(x, np.real(fn),  linestyle='-.', color="r", label='Re(fn)')
plt.plot(x, np.imag(fn),  linestyle=':', color="g", label='Im(fn)')
#plt.plot(x, fa, linestyle='solid', color="y", label='fa')
plt.xlabel('xi')
plt.ylabel('F')
plt.legend()
plt.axis([-50, 50, 0, 0.1])
plt.show()
