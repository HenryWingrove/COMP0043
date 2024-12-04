clear all
close all
format long

% Set the parameters
ngrid = 2^(10); % number of grid points
% distr = 6; % kernel: 1 = normal, ..., 8 = stable; see parameters.m
% ndates = 252; % monitoring dates
% Euler = 1; % 0 = never use Euler acceleration
% T = 1; % maturity
% dt = T/ndates; % time step
% rf = 0.05; % risk free interest rate
% q = 0.02; % dividend rate
% S_0 = 100; % spot price
% K = 90; % strike price
distr = 1; % kernel: 1 = normal, ..., 8 = stable; see parameters.m
ndates = 25; % monitoring dates
Euler = 1; % 0 = never use Euler acceleration
T = 0.5; % maturity
dt = T/ndates; % time step
rf = 0.1; % risk free interest rate
q = 0; % dividend rate

% Compute the parameters for the grid bounds and the payoff
param = parameters(distr,T,dt,rf,q);
[low,up] = bounds(param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Spitzer with sinc')
% BASATO SU Green, Fusai, Abrahams (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0; % shift parameter = damping factor

b = up; % upper bound of the support
trunc_L = exp(-b); % lower truncation
trunc_U = exp(b); % upper truncation

[x,h,xi,H] = kernel(ngrid,-b,b,param,alpha,1,1);

% sol = solve_allq_FHT_S_LB_dist(H,b,ndates,Euler,@(x)ifht(x,'sinc',ngrid),-1,1);
% distribuzione1 = sol;
% sol = solve_allq_FHT_S_LB_dist(H,b,ndates,Euler,@(x)ifht(x,'sinc',ngrid),-1,2);
% distribuzione2 = sol;
sol = solve_allq_FHT_S_LB_dist(H,b,ndates,Euler,@(x)ifht(x,'sgn0',0),-1,1);
distribuzione1 = sol;
sol = solve_allq_FHT_S_LB_dist(H,b,ndates,Euler,@(x)ifht(x,'sgn0',0),-1,2);
distribuzione2 = sol;

figure
plot(x,distribuzione1);
hold on
plot(x,distribuzione2,'r:');
figure
plot(x,distribuzione1-distribuzione2);

