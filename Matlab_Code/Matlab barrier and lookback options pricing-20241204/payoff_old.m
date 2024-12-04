% Compute the scale, the payoff and its Fourier transform
function [S,g,G] = payoff(x,xi,alpha,K,L,U,C,call)

% Scale
S = C*exp(x);

% Payoff
%ind_fun = (S>=L).*(S<=U);
%if call == 1 % see also GFA, Eq. (3.5)
%    g = exp(alpha*(x+log(C/K))).*max(S-K,0).*ind_fun;
%else % put
%    g = exp(alpha*(x+log(C/K))).*max(K-S,0).*ind_fun;
%end
g = exp(alpha*(x+log(C/K))).*max(call*(S-K),0).*(S>=L).*(S<=U);

% Analytical Fourier transform of the payoff
if call == 1 % call
    b = log(U/K);
else % put
    b = log(L/K);
end
% L. Feng, V. Linetsky, Math. Finance 18, 337-384 (2008), Eq. (4.7)
G = K*((1-exp(b*(alpha+1i*xi)))./(alpha+1i*xi) ...
   - (1-exp(b*(1+alpha+1i*xi)))./(1+alpha+1i*xi));
% L. Feng, V. Linetsky, Math. Finance 18, 337-384 (2008), Eq. (4.6)
%G = -K./((-1i*alpha+xi).*(-1i*(1+alpha)+xi));
% GFA (2010), Eq. (3.7b)
%b = max(log(L/C),log(K/C));
%G = C*(exp(log(K/C)+b*(alpha+1i*xi))./(alpha+1i*xi) ...
%    - exp(b*(1+alpha+1i*xi))./(1+alpha+1i*xi));
if (alpha == 0) % eliminable discontinuity; otherwise 0/0 = NaN
    G(floor(end/2)+1) = K*(exp(b)-b-1);
end
% Adapt Eq. (4.7), which is written for S = K*exp(x), to S = C*exp(x) using
% the shift theorem g(x+a) <-> G(xi)*exp(-1i*a*xi) with a = log(C/K)
if C ~= K
    G = G.*exp(log(K/C)*(alpha+1i*xi));  
end
% Numerico
% G = fftshift(ifft(ifftshift(g)))*(x(2)-x(1))*length(x);

% % Comparison of the numerical and analytical Fourier transform of the payoff
% 
% % Normal space
% gn = fftshift(fft(ifftshift(G)))./((x(2)-x(1))*length(x));
% figure(1), clf
% plot(x,g,'ro',x,real(gn),'gs')
% xlabel('x')
% ylabel('Re(g)')
% legend('analytical','numerical')
% figure(2), clf
% plot(x,zeros(size(x)),'ro',x,imag(gn),'gs')
% xlabel('x')
% ylabel('Im(g)')
% legend('analytical','numerical')
% figure(3), clf
% plot(x,g./real(gn),'gs')
% xlabel('x')
% ylabel('ga/Re(gn)')
% xlim([0 log(U/C)])
% figure(4), clf
% plot(x,g,'r',x,real(gn),':r',x,imag(gn),':g')
% legend('ga','Re(gn)','Im(gn)')
%
% % Reciprocal space
% Gn = fftshift(ifft(ifftshift(g)))*(x(2)-x(1))*length(x);
% figure(5), clf
% plot(xi,real(G),'ro',xi,real(Gn),'gs')
% xlabel('xi')
% ylabel('Re(G)')
% legend('analytical','numerical')
% xlim([-50 50])
% figure(6), clf
% plot(xi,imag(G),'ro',xi,imag(Gn),'gs')
% xlabel('xi')
% ylabel('Im(G)')
% legend('analytical','numerical')
% xlim([-50 50])
% figure(7), clf
% plot(xi,real(G)./real(Gn),'ro',xi,imag(G)./imag(Gn),'gs')
% xlabel('xi')
% legend('Re(Ga)/Re(Gn)','Im(Ga)/Im(Gn)')
% xlim([-50 50])
% figure(8), clf
% plot(xi,real(G),'r',xi,real(Gn),':r',xi,imag(G),'g',xi,imag(Gn),':g')
% legend('Re(Ga)','Re(Gn)','Im(Ga)','Im(Gn)')
% xlim([-50 50])