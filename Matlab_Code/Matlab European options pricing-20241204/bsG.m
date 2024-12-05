% Parameters
xwidth = 8;        % Width of the support in real space
ngrid_real = 2^9;  % Finer grid size in log-price space (512 points)
dx_real = xwidth / ngrid_real;  % Step size in real space
x_real = linspace(-xwidth/2, xwidth/2, ngrid_real);  % Log-price grid

% PDF of arithmetic Brownian motion at maturity
muABM = r - q - 0.5 * sigma^2; % Drift coefficient
pdf_XT = @(x) exp(-(x - muABM * T).^2 / (2 * sigma^2 * T)) / sqrt(2 * pi * sigma^2 * T);

% Undamped payoff function
g_call = max(S0 * exp(x_real) - K, 0);  % Call
g_put = max(K - S0 * exp(x_real), 0);   % Put

% Numerical quadrature for the second line
V_call_quad = exp(-r * T) * trapz(x_real, g_call .* pdf_XT(x_real));
V_put_quad = exp(-r * T) * trapz(x_real, g_put .* pdf_XT(x_real));

% Display results
fprintf('Numerical Quadrature:\n');
fprintf('Call Price: %.10f\n', V_call_quad);
fprintf('Put Price: %.10f\n', V_put_quad);