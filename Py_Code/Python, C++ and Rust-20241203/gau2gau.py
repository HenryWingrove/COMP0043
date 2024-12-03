# terrible style - NO GLOBAL VARIABLES!

import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

N = 1024
Dx = 0.1
x = Dx*np.arange(-N/2, N/2)
L = N*Dx
Dxi = 2*np.pi/L
xi = Dxi*np.arange(-N/2, N/2)
W = N*Dxi

a = 2
fa = np.sqrt(a/np.pi) * np.exp(-a*np.power(x, 2))
Fa = np.exp(-np.power(xi, 2)/(4*a))

Fn = fftshift(ifft(ifftshift(fa)))*L
fn = fftshift(fft(ifftshift(Fa)))/L
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
