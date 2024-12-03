# terrible style - NO GLOBAL VARIABLES!

import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

N = 2024
Dx = 0.01
Lx = N*Dx
Dxi = 2*np.pi/Lx
Lxi = N*Dxi
x = Dx*np.arange(-N/2, N/2)
xi = Dxi*np.arange(-N/2, N/2)

a = 1
fa = a/2*np.exp(-a*abs(x))
Fa = np.divide(pow(a, 2), pow(a, 2)+np.power(xi, 2))

Fn = fftshift(ifft(ifftshift(fa)))*Lx
fn = fftshift(fft(ifftshift(Fa)))/Lx
Fn1 = fftshift(fft(ifftshift(fa)))*Dx
fn1 = fftshift(ifft(ifftshift(fa)))/Dx

plt.figure(1)
plt.plot(x, np.real(fn),  linestyle='-.', color="r", label='Re(fn)')
plt.plot(x, np.imag(fn),  linestyle=':', color="g", label='Im(fn)')
plt.plot(x, fa, linestyle='solid', color="y", label='fa')
plt.xlabel('xi')
plt.ylabel('F')
plt.legend()
plt.axis([-5, 5, 0, 1])
plt.show()
