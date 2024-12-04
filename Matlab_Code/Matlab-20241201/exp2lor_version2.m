%% Check numerically the Fourier pair Laplace <-> Lorentzian

N = 1024;
Dx = 0.1;
x = Dx*(-N/2:N/2-1);
Lx = N*Dx;
% Nyquist relation: Dx*Dnu = 1/N
Dnu = 1/Lx;
nu = 2*pi*Dnu*(-N/2:N/2-1);
Nu = N*Dnu; % = 1/Dt;

a = 1;
fa = a/2*exp(-a*abs(x)); % Laplace (or double-sided exponential)
Fa = a^2./(a^2+(nu).^2); % Lorentzian (or Cauchy)

Fn  = fftshift(ifft(ifftshift(fa)))*Lx;
fn  = fftshift( fft(ifftshift(Fa)))/Lx;
Fn1 = fftshift( fft(ifftshift(fa)))*Dx;
fn1 = fftshift(ifft(ifftshift(Fa)))/Dx;

close all
figure(5), clf, hold on
plot(x,real(fn),'r')
plot(x,imag(fn),'g')
plot(x,fa,'y:')
axis([-10 10 0 1])
xlabel('t')
ylabel('f')
legend('Re(fn)','Im(fn)','fa')

figure(6), clf, hold on
plot(nu,real(Fn),'b')
plot(nu,imag(Fn),'m')
plot(nu,Fa,'c:')
axis([-10/pi 10/pi 0 1])
xlabel('\nu')
ylabel('F')
legend('Re(Fn)','Im(Fn)','Fa')

figure(7), clf, hold on
plot(x,real(fn1),'r')
plot(x,imag(fn1),'g')
plot(x,fa,'y:')
axis([-10 10 0 1])
xlabel('t')
ylabel('f')
legend('Re(fn1)','Im(fn1)','fa')

figure(8), clf, hold on
plot(nu,real(Fn1),'b')
plot(nu,imag(Fn1),'m')
plot(nu,Fa,'c:')
axis([-10/pi 10/pi 0 1])
xlabel('\omega')
ylabel('F')
legend('Re(Fn1)','Im(Fn1)','Fa')

