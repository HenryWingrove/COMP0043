% Price a lookback option using Spitzer identity
clear all
close all
format long

for i = 8:14
% Set the parameters
ngrid = 2^i; % number of grid points
distr = 1; % kernel: 1 = normal, ..., 8 = stable; see parameters.m
ndates = 50; % monitoring dates
Euler = 1; % 0 = never use Euler acceleration
T = 0.5; % maturity
dt = T/ndates; % time step
rf = 0.1; % risk free interest rate
q = 0; % dividend rate
S_0 = 1; % spot price
K = 1; % strike price
call = 1; % 1 = call, -1 = put

% Compute the parameters for the grid bounds and the payoff
param = parameters(distr,T,dt,rf,q);
[low,up] = bounds(param);
%b = max(abs(low),up); % upper bound of the support
b = up; % upper bound of the support
trunc_L = S_0*exp(-b); % lower truncation
trunc_U = S_0*exp(b); % upper truncation

alpha = 0; % shift parameter = damping factor
[x,h,xi,H] = kernel(ngrid,-b,b,param,alpha,1,1);
C=S_0;
%[S,g,G] = payoff(x,xi,alpha,K,trunc_L,trunc_U,C,call);
[S,g,G] = payoff_old(x,xi,alpha,K,trunc_L,trunc_U,C,call);

%disp('Spitzer with sinc')
% Basato su Green, Fusai, Abrahams (2010), Eq. (3.21)
%G = -K*exp(1i*xi*log(K/S_0))./(xi.^2-1i*xi);
tic
sol = solve_allq_FHT_S_LB(H,G,b,ndates,Euler,@(x)ifht(x,'sinc',ngrid),call,x,C,S_0);
price = S_0*sol;
CPU_Time = toc;

fprintf('%2d %14.10f %14.10f\n',i,price,CPU_Time)
end