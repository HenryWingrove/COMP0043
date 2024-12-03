%% Check numerically the Fourier pair Laplace <-> Lorentzian

% Grids in real and Fourier space
N = 4096;                % Increase grid size for better resolution
Dx = 0.01;               % Grid step in real space
Lx = N * Dx;             % Upper truncation limit in real space
Dxi = 2 * pi / Lx;       % Grid step in Fourier space
x = Dx * (-N/2:N/2-1);   % Grid in real space
xi = Dxi * (-N/2:N/2-1); % Grid in Fourier space

% Analytical expressions
lambda = 3;                              % Activity or rate parameter
fa = exp(-lambda * abs(x));              % Laplace (real space)
Fa = 2 * lambda ./ (xi.^2 + lambda^2);   % Lorentzian (Fourier space)

% Numerical Fourier Transform: Properly Scaled
Fn = fftshift(fft(ifftshift(fa))) * Dx;   % Fourier Transform of fa -> Fourier space
fn = fftshift(ifft(ifftshift(Fa))) * Dxi; % Inverse Fourier Transform of Fa -> Real space

% Check scaling
disp(['Sum of fa: ', num2str(sum(fa) * Dx)]);
disp(['Sum of real(fn): ', num2str(sum(real(fn)) * Dx)]);

% Plot 1: Real-space function (fa) and its numerical inverse FFT
figure(1), clf, hold on
plot(x, fa, 'k:', 'LineWidth', 2, 'DisplayName', 'Analytical f(x)')
plot(x, real(fn), 'r', 'DisplayName', 'Re IFFT(F)')
plot(x, imag(fn), 'g', 'DisplayName', 'Im IFFT(F)')
axis([-4 4 -0.2 1.2])  % Adjust scale
xlabel('x'), ylabel('f(x)')
legend('show')
title('Laplace or double exponential function')

% Plot 2: Fourier-space function (Fa) and its numerical FFT
figure(2), clf, hold on
plot(xi, Fa, 'c:', 'LineWidth', 2, 'DisplayName', 'Analytical F(\xi)')
plot(xi, real(Fn), 'b', 'DisplayName', 'Re FFT(f)')
plot(xi, imag(Fn), 'm', 'DisplayName', 'Im FFT(f)')
axis([-20 20 -0.2 1.2])  % Adjust scale
xlabel('\xi'), ylabel('F(\xi)')
legend('show')
title('Fourier Transform: Analytical vs Numerical')

% Plot 3: Real-space function (fa) and its numerical inverse FFT
figure(3), clf, hold on
plot(x, fa, 'k:', 'LineWidth', 2, 'DisplayName', 'Analytical f(x)')
plot(x, real(fn), 'r', 'DisplayName', 'Re IFFT(F)')
plot(x, imag(fn), 'g', 'DisplayName', 'Im IFFT(F)')
axis([-4 4 -0.2 1.2])  % Adjust scale
xlabel('x'), ylabel('f(x)')
legend('show')
title('Laplace function with corrected Inverse FFT')

% Plot 4: Fourier-space function (Fa) and its numerical scaled FFT
figure(4), clf, hold on
plot(xi, Fa, 'c:', 'LineWidth', 2, 'DisplayName', 'Analytical F(\xi)')
plot(xi, real(Fn1), 'b', 'DisplayName', 'Re FFT(f) (scaled)')
plot(xi, imag(Fn1), 'm', 'DisplayName', 'Im FFT(f) (scaled)')
axis([-20 20 -0.2 1.2])  % Adjust scale
xlabel('\xi'), ylabel('F(\xi)')
legend('show')
title('Fourier Transform with scaled FFT')