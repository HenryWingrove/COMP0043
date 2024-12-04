% Price a single or double barrier option using fast Fourier transform methods

%clc
clear all
format long

% Set the parameters
ngrid = 2^14; % number of grid points
distr = 2; % kernel: 1 = normal, ..., 8 = stable; see parameters.m
ndates = 52; % monitoring dates
Euler = 1; % 0 = do not use Euler acceleration
T = 1; % maturity
dt = T/ndates; % time step
rf = 0.05; % risk-free interest rate
q = 0.02; % dividend rate
S_0 = 1; % spot price
K = 0.6; % strike price
call = 0; % 1 for call, other number for put

% fattore di conversione 
C=S_0;

% Compute the parameters for the grid bounds and for the payoff
param = parameters(distr,T,dt,rf,q);
[low,up] = bounds(param);
up=max(abs(low),up);
trunc = S_0*exp(up); % truncation
barrier_L = 0.8; % lower barrier
barrier_U = trunc; % down-and-out
b = log(trunc/C); % upper bound of the support

%alpha = -14; % shift parameter
alpha = set_alpha(param,call,log(barrier_L/C),log(barrier_U/C),b)

% Print header
fprintf('\nndates = %i, ngrid = %d\n',ndates,ngrid)
fprintf('Method           Price         CPU_t/s      Average_iterations\n')

% Compute the grid and the discounted kernel for Fourier-based methods
[x,h,xi,H] = kernel(ngrid,-b,b,param,alpha,1,0);

% Compute the scale, payoff and its FT for Fourier-based methods
[S,g,G] = payoff(x,xi,alpha,K,barrier_L,barrier_U,C,call);
%[S,g,G] = payoff_old(x,xi,alpha,K,barrier_L,barrier_U,C,call);

% ripeto a vuoto per misurare il tempo
[sol,miter] = solve_allq_FHT_DB(H,G,log(barrier_L/C),log(barrier_U/C),b,ndates,Euler,@(x)ifht(x,'sinc',ngrid));
sol = sol.*exp(-alpha*x);
price = 100*interp1(S,sol,S_0,'spline');
[sol,miter] = solve_allq_FHT_DB(H,G,log(barrier_L/C),log(barrier_U/C),b,ndates,Euler,@(x)ifht(x,'sinc',ngrid));
sol = sol.*exp(-alpha*x);
price = 100*interp1(S,sol,S_0,'spline');
%%%%%%%%%%%%%%5

tic
[sol,miter] = solve_allq_FHT_DB(H,G,log(barrier_L/C),log(barrier_U/C),b,ndates,Euler,@(x)ifht(x,'sinc',ngrid));
sol = sol.*exp(-alpha*x);
price = 100*interp1(S,sol,S_0,'spline');
fprintf('Z-WH         %10.12f    %9f    %f\n',price,toc,miter)

tic
sol = feng_linetsky(H,G,log(barrier_L/C),log(barrier_U/C),b,ndates,@(x)ifht(x,'sinc',ngrid));
sol = sol.*exp(-alpha*x);
price = 100*interp1(S,sol,S_0,'spline');
fprintf('Feng-Linetsky    %10.12f    %9f\n',price,toc)

tic
% coniugo H perchè forward
sol = solve_allq_FHT_S(conj(H),G,log(barrier_L/C),b,ndates,Euler,@(x)ifht(x,'sinc',ngrid));
price = sol(end/2+1)*100;
fprintf('Z-Spitzer    %10.12f    %9f    %f\n',price,toc,miter)
