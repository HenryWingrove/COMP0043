%% Compute the PDF of non-central chi-squared distribution via inverse Fourier transform

% Original script parameters
N = 1024;
Dx = 0.1;
x = Dx*(-N/2:N/2-1);
Lx = N*Dx;
% Nyquist relation: Dx*Dnu = 1/N
Dnu = 1/Lx;
nu = 2*pi*Dnu*(-N/2:N/2-1);
Nu = N*Dnu; % = 1/Dt;

% Parameters of the non-central chi-squared distribution
k = 5;        % Degrees of freedom
lambda = 5;   % Non-centrality parameter

% Analytical PDF of the non-central chi-squared distribution (fa)
fa = zeros(size(x));
for i = 1:length(x)
    xi = x(i);
    if xi >= 0
        fa(i) = 0.5 * exp(-(xi + lambda)/2) * (xi / lambda)^(k/4 - 0.5) * ...
                 besseli(k/2 - 1, sqrt(lambda * xi));
    else
        fa(i) = 0;
    end
end

% Characteristic function of the non-central chi-squared distribution (Fa)
Fa = exp(1j * lambda * nu) .* exp(-1j * nu / 2) ./ (1 - 2j * nu).^(k/2);

% Compute Fourier transforms using the original script's conventions
Fn  = fftshift(ifft(ifftshift(fa)))*Lx;   % Fourier transform of fa
fn  = fftshift( fft(ifftshift(Fa)))/Lx;   % Inverse Fourier transform of Fa

% Plotting results (following the original script's plotting conventions)
close all

% Plot fa and fn to compare the analytical PDF and numerical inverse Fourier transform
figure(1), clf, hold on
plot(x, real(fn), 'r', 'LineWidth', 1.5)
plot(x, imag(fn), 'g', 'LineWidth', 1.5)
plot(x, fa, 'b--', 'LineWidth', 1.5)
axis([-5 30 0 max(fa)*1.1])
xlabel('x')
ylabel('PDF')
legend('Re(fn)', 'Im(fn)', 'Analytical fa')
title('Non-central Chi-squared PDF via Inverse Fourier Transform')
grid on

% Plot Fa and Fn to compare the characteristic functions
figure(2), clf, hold on
plot(nu, real(Fn), 'b', 'LineWidth', 1.5)
plot(nu, imag(Fn), 'm', 'LineWidth', 1.5)
plot(nu, real(Fa), 'c--', 'LineWidth', 1.5)
axis([-10 10 -0.5 1])
xlabel('\nu')
ylabel('Characteristic Function')
legend('Re(Fn)', 'Im(Fn)', 'Re(Fa)')
title('Characteristic Function of Non-central Chi-squared Distribution')
grid on
