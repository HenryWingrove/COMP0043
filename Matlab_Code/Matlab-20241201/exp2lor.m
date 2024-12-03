% Grids in real and Fourier space
% Linked by the Nyquist relation Dx*Dxi = 2*pi/N
N = 2048;               % Grid size
Dx = 0.01;              % Grid step in real space
Lx = N * Dx;            % Upper truncation limit in real space
Dxi = 2 * pi / Lx;      % Grid step in Fourier space
Lxi = N * Dxi;          % Upper truncation limit in Fourier space
x = Dx * (-N/2:N/2-1);  % Grid in real space
xi = Dxi * (-N/2:N/2-1); % Grid in Fourier space

% Analytical expressions
lambda = 3;             % Activity or rate parameter
fa = exp(-lambda * abs(x));          % Laplace or double exponential function
Fa = 2 * lambda ./ (xi.^2 + lambda^2); % Lorentz or Cauchy function

% Numerical approximations
Fn = fftshift(ifft(ifftshift(fa))) * Lx; % Fourier transform of fa
fn = fftshift(fft(ifftshift(Fa))) / Lx; % Inverse Fourier transform of Fa
Fn1 = fftshift(fft(ifftshift(fa))) * Dx; % Fourier transform with Dx scaling
fn1 = fftshift(ifft(ifftshift(Fa))) / Dx; % Inverse Fourier with Dx scaling

% Close all figures to start fresh
close all

% Figure 1: Laplace function and numerical IFFT
figure(1), clf, hold on
plot(x, fa, 'r', 'DisplayName', 'Analytic f(x)')
plot(x, real(fn), 'g.', 'DisplayName', 'Re IFFT(F)')
plot(x, imag(fn), 'b.', 'DisplayName', 'Im IFFT(F)')
axis([-4 4 0 1.2])
xlabel('x')
ylabel('f')
legend
title('Laplace or double exponential function')

% Additional Figures

% Figure 2: Fourier transform visualization
figure(2), clf, hold on
plot(xi, Fa, 'r', 'DisplayName', 'Analytical Fa')
plot(xi, real(Fn), 'g.', 'DisplayName', 'Re FFT(f)')
plot(xi, imag(Fn), 'b.', 'DisplayName', 'Im FFT(f)')
xlabel('xi')
ylabel('F')
legend
title('Fourier Transform of Laplace Function')
axis([-50 50 -0.5 1.5])

% Figure 3: Amplitude comparison
figure(3), clf, hold on
plot(xi, Fa, 'r', 'DisplayName', 'Analytical Fa')
plot(xi, abs(Fn), 'g--', 'DisplayName', '|FFT(f)|')
xlabel('xi')
ylabel('Amplitude')
legend
title('Analytical vs Numerical Amplitudes')
axis([-50 50 0 1.5])

% Figure 4: Reconstruction of original function
figure(4), clf, hold on
plot(x, fa, 'r', 'DisplayName', 'Analytical f(x)')
plot(x, real(fn1), 'g.', 'DisplayName', 'Re IFFT(Fn1)')
plot(x, imag(fn1), 'b.', 'DisplayName', 'Im IFFT(Fn1)')
xlabel('x')
ylabel('f')
legend
title('Reconstruction from Fourier Space')
axis([-4 4 0 1.2])