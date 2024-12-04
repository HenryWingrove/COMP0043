% Compute the scale, the payoff and its Fourier transform
function [S,g,G] = payoff(x,xi,alpha,K,L,U,C,call)

% Scale
S = C*exp(x);

% Payoff 
if call == 1 % see e.g. Green, Fusai, Abrahams 2010, Eq. (3.24)
    g = exp(alpha*x).*max(S-K,0).*(S>=L).*(S<=U);
else % put
    g = exp(alpha*x).*max(K-S,0).*(S>=L).*(S<=U);
end

% Analytical Fourier transform of the payoff

% Integration bounds
if call == 1 % call
    b1 = log(max(L,K)/C);
    b2 = log(U/C);
else % put
    b1 = log(min(U,K)/C);
    b2 = log(L/C);
end

xi2 = alpha+1i*xi;
k = log(K/C);

% GFA 2010 Eq. (3.26) with extension to put option
G = C*((exp(b2*(1+xi2))-exp(b1*(1+xi2)))./(1+xi2) ...
    - (exp(k+b2*xi2)-exp(k+b1*xi2))./xi2);

% Eliminable discontinuities for xi = 0; otherwise 0/0 = NaN
if (alpha == 0)
    G(floor(end/2)+1) = C*(exp(b2)-exp(b1)-exp(k)*(b2-b1));
elseif (alpha == -1)
    G(floor(end/2)+1) = C*(b2-b1+exp(k-b2)-exp(k-b1));
end

% Plot to compare the analytical and numerical payoffs
gn = fftshift(fft(ifftshift(G)))./((x(2)-x(1))*length(x));
figure(1), clf
plot(x,g,'r',x,real(gn),'g')
xlabel('x')
ylabel('Re(g)')
legend('analytical','numerical')
title('Payoff function')