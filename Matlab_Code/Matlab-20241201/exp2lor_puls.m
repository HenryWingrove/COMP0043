%% Check numerically the Fourier pair Laplacian <-> Lorentzian

% Grids in real and Fourier space, linked by the Nyquist relation
% Dx*Dxi = 2*pi/N
N = 2048; % grid size
Dx = 0.01; % grid step in real space
Lx = N*Dx; % upper truncation limit in real space
Dxi = 2*pi/Lx; % grid step in Fourier space
Lxi = N*Dxi; % upper truncation limit in Fourier space
x = Dx*(-N/2:N/2-1); % grid in real space
xi = Dxi*(-N/2:N/2-1); % grid in Fourier space

% Analytical expressions
fa = 1/(2*pi)*exp(-pi*abs(x)); % Laplace or bilateral exponential function
Fa = 1./(pi^2+xi.^2); % Lorentz or Cauchy function

% Numerical approximations
Fn = fftshift(ifft(ifftshift(fa)))*Lx;
fn = fftshift(fft(ifftshift(Fa)))/Lx;

close all
figure(1), clf, hold on
xlabel('\xi')
ylabel('F')
axis([-8 8 0 0.16])
plot(xi,Fa,'r')
plot(xi,real(Fn),'ob')
plot(xi,imag(Fn),':g')
legend('Fa','Re(Fn)','Im(Fn)')
title('Lorentz or Cauchy function')

figure(2), clf, hold on
xlabel('x')
ylabel('f')
axis([-8 8 0 0.16])
plot(x,fa,'r')
plot(x,real(fn),'ob')
plot(x,imag(fn),':g')
legend('fa','Re(fn)','Im(fn)')
title('Laplace or double exponential function')