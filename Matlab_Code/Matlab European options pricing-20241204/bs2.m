%% Compute the Black-Scholes-Merton price of European options
% The model for the underlying is geometric Brownian motion
% dS = mu*S*dt + sigma*S*dW

clear variables % clear all variables from the workspace

% Contract parameters
T = 1; % maturity
K = 1.1; % strike price

% Market parameters
S0 = 1; % initial stock price
r = 0.05; % risk-free interest rate
q = 0.02; % dividend rate

% Model parameter
sigma = 0.4; % volatility

% Direct integration parameters
xwidth = 8; % width of the support in real space
ngridx = 2^16; % number of grid points

% Fourier parameters
ngrid = 2^6; % number of grid points
alphac = -6; % damping parameter for call
alphap = 6; % damping parameter for put

% Monte Carlo parameters; paths = nblocks*npaths
nblocks = 10000; % number of blocks
npaths = 2000; % number of paths per block

% Controls
figures = 0;

%% Analytical solution

% Black-Scholes, Matlab's Financial Toolbox
tic
[VcaM,VpaM] = blsprice(S0,K,r,T,sigma,q);
cputime_aM = toc;

% Black-Scholes, own
tic
muABM = r-q-0.5*sigma^2; % drift coefficient of the arithmetic Brownian motion
d2 = (log(S0/K)+muABM*T)/(sigma*sqrt(T));
d1 = d2 + sigma*sqrt(T);
Vcao = S0*exp(-q*T)*cdf('Normal',d1,0,1) - K*exp(-r*T)*cdf('Normal',d2,0,1);
Vpao = K*exp(-r*T)*cdf('Normal',-d2,0,1) - S0*exp(-q*T)*cdf('Normal',-d1,0,1);
% Put-call parity: Vp + S0*exp(-q*T) = Vc + K*exp(-r*T) 
cputime_ao = toc;

% Discounted expected payoff 
tic
Nx = ngridx/2;
dx = xwidth/ngridx; % grid step of the log-price
x = dx*(-Nx:Nx-1); % grid of the log-price
f = 1./(sqrt(2*pi*T)*sigma).*exp(-(x-muABM*T).^2./(2*sigma^2*T)); % PDF at maturity for ABM
Vcae = exp(-r*T)*sum(max(S0*exp(x)-K,0).*f)*dx; % call
Vpae = exp(-r*T)*sum(max(K-S0*exp(x),0).*f)*dx; % put
cputime_ae = toc;

gc = max(S0*exp(x)-K,0).*f;
max(abs(gc(2:end)-gc(1:end-1)))
gp = max(K-S0*exp(x),0).*f;
max(abs(gp(2:end)-gc(1:end-1)))

% Print the results
fprintf('%-20s%15s%15s%15s\n','','Call','Put','CPU_time/s')
fprintf('%-20s%15.10f%15.10f%15.10f\n','Black-Scholes Matlab',VcaM,VpaM,cputime_aM)
fprintf('%-20s%15.10f%15.10f%15.10f\n','Black-Scholes own',Vcao,Vpao,cputime_ao)
fprintf('%-20s%15.10f%15.10f%15.10f\n','Disc.expected payoff',Vcae,Vpae,cputime_ae)

if figures

    % Plot the analytical solution
    [St,t] = meshgrid(0:.05:2,0:0.025:T);
    d2 = (log(St/K)+muABM*(T-t))./(sigma*sqrt(T-t));
    d1 = d2 + sigma*sqrt(T-t);

    close all
    figure(1)
    Vc = St.*exp(-q*(T-t)).*cdf('Normal',d1,0,1) - K*exp(-r*(T-t)).*cdf('Normal',d2,0,1);
    Vc(end,:) = max(St(end,:)-K,0);
    mesh(St,t,Vc)
    xlabel('S')
    ylabel('t')
    zlabel('V')
    title('Call')
    view(-30,24)
    print('-dpng','bsc.png')
    
    figure(2)
    Vp = K*exp(-r*(T-t)).*cdf('Normal',-d2,0,1) - St.*exp(-q*(T-t)).*cdf('Normal',-d1,0,1);
    Vp(end,:) = max(K-St(end,:),0);
    mesh(St,t,Vp)
    xlabel('S')
    ylabel('t')
    zlabel('V')
    title('Put')
    view(30,24)
    print('-dpng','bsp.png')
    
    % Plot the analytical solution as a function of the log price
    k = log(K/S0); % log-strike
    [xt,t] = meshgrid(-1:.05:1,0:0.025:T);
    d2 = (xt-k+muABM*(T-t))./(sigma*sqrt(T-t));
    d1 = d2 + sigma*sqrt(T-t);
    
    figure(3)
    Vc = S0*(exp(xt-q*(T-t)).*cdf('Normal',d1,0,1)-exp(k-r*(T-t)).*cdf('Normal',d2,0,1));
    Vc(end,:) = S0*max(exp(xt(end,:))-exp(k),0);
    mesh(xt,t,Vc)
    xlabel('x')
    ylabel('t')
    zlabel('V')
    title('Call')
    view(-30,24)
    print('-dpng','bscx.png')
    
    figure(4)
    Vp = S0*(exp(k-r*(T-t)).*cdf('Normal',-d2,0,1)-exp(xt-q*(T-t)).*cdf('Normal',-d1,0,1));
    Vp(end,:) = S0*max(exp(k)-exp(xt(end,:)),0);
    mesh(xt,t,Vp)
    xlabel('x')
    ylabel('t')
    zlabel('V')
    title('Put')
    view(30,24)
    print('-dpng','bspx.png')

end

%% Fourier transform method

% Grid in Fourier space
tic
N = ngrid/2;
dx = xwidth/ngrid; % grid step of the log-price
x = dx*(-N:N-1); % grid of the log-price
b = xwidth/2; % upper bound of the log-price
dxi = pi/b; % Nyquist relation: grid step in Fourier space
xi = dxi*(-N:N-1); % grid in Fourier space

% These functions provide the characteristic functions of 8 Levy processes
% param = parameters(1,T,T,r,q); % set the parameters editing parameters.m
% [x,fc,xi,Psic] = kernel(ngrid,-b,b,param,alphac,0,1); % call
% [x,fp,xi,Psip] = kernel(ngrid,-b,b,param,alphap,0,1); % put

% Damped payoff and its Fourier transform
U = S0*exp(b); % upper bound of the spot price
L = S0*exp(-b); % lower bound of the spot price
[~,gc,Gc] = payoff(x,xi,alphac,K,L,U,S0,1,figures); % call
[S,gp,Gp] = payoff(x,xi,alphap,K,L,U,S0,-1,figures); % put

% Characteristic function at maturity for arithmetic Brownian motion
psi = @(xi) 1i*muABM*xi-0.5*(sigma*xi).^2; % characteristic exponent
Psic = exp(psi(xi+1i*alphac)*T); % shifted characteristic function for a call
Psip = exp(psi(xi+1i*alphap)*T); % shifted characteristic function for a put

% Discounted expected payoff computed with the Plancherel theorem
% Because of symmetry, it is sufficient to integrate on the half grid
%c = exp(-r*T).*real(fftshift(fft(ifftshift(Gc.*conj(Psic)))))/xwidth; % call
%VcF = interp1(S,c,S0,'spline');
VcF = exp(-r*T)/pi*trapz(real(Gc(N+1:end).*conj(Psic(N+1:end))))*dxi; % call
%p = exp(-r*T).*real(fftshift(fft(ifftshift(Gp.*conj(Psip)))))/xwidth; % put
%VpF = interp1(S,p,S0,'spline');
VpF = exp(-r*T)/pi*trapz(real(Gp(N+1:end).*conj(Psip(N+1:end))))*dxi; % put
cputime_F = toc;
fprintf('%-20s%15.10f%15.10f%15.10f\n','Fourier',VcF,VpF,cputime_F)

hc = real(Gc(N+1:end).*conj(Psic(N+1:end)));
max(abs(hc(2:end)-hc(1:end-1)))
hp = real(Gp(N+1:end).*conj(Psip(N+1:end)));
max(abs(hp(2:end)-hp(1:end-1)))

% figures_ft(S,x,xi,f,Psic,gc,Gc) % call

%% Monte Carlo

tic;
VcMCb = zeros(nblocks,1);
VpMCb = zeros(nblocks,1);
for i = 1:nblocks

    % Arithmetic Brownian motion X(T) = log(S(T)/S(0)) at time T
    X = muABM*T + sigma*randn(1,npaths)*sqrt(T);

    % Transform to geometric Brownian motion S(T) at time T
    S = S0*exp(X);

    % Discounted expected payoff
    VcMCb(i) = exp(-r*T)*mean(max(S-K,0));
    VpMCb(i) = exp(-r*T)*mean(max(K-S,0));

end
VcMC = mean(VcMCb);
VpMC = mean(VpMCb);
scMC = sqrt(var(VcMCb)/nblocks);
spMC = sqrt(var(VpMCb)/nblocks);
cputime_MC = toc;
fprintf('%-20s%15.10f%15.10f%15.10f\n','Monte Carlo',VcMC,VpMC,cputime_MC)
fprintf('%-20s%15.10f%15.10f\n','Monte Carlo stdev',scMC,spMC)
