% Price a Lookback Option using several quadrature or fast Fourier transform methods
clear all
close all
format long

fprintf('ZT+QuadTrap ------ ZT+Quad4th ------- ZT+FHT ----------- Spitzer\n');
fprintf('Prezzo - Tempo --- Prezzo - Tempo --- Prezzo - Tempo --- Prezzo - Tempo\n');

for i=1:6
% Set the parameters
ngrid = 2^(10+i); % number of grid points
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
ndates = 50; % monitoring dates
Euler = 1; % 0 = never use Euler acceleration
T = 0.5; % maturity
dt = T/ndates; % time step
rf = 0.1; % risk free interest rate
q = 0; % dividend rate
S_0 = 85; % spot price
K = 90; % strike price

alpha=0;

% Compute the parameters for the grid bounds and the payoff
param = parameters(distr,T,dt,rf,q);
[low,up] = bounds(param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('ZT Algorithms') 
% BASATI SU RICORSIONE Fusai, Marazzina, Marena, Ng (2012)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx = log(S_0);     lk = log(K);
l = 0;             u = lx+up;

b = u-l; % upper bound of the support
%kernel_LB as kernel without the discounting factor
[x,h,w,H] = kernel(ngrid,-b,b,param,alpha,0,0);
% Compute scale, payoff and its FT 
N = ngrid/2;
g = exp(alpha*x).*[zeros(1,N+1) ones(1,N-1)]; % payoff
G = fftshift(ifft(ifftshift(g)))*2*b;

x = x(ngrid/2+1:end);
Ris(1,1:8)=0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LSflag = 2; % 0 = Gaussian elimination, 1 = GMRES, 2 = GMRES with preconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('- Quadrature with trapz (2nd order method)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
sol = solve_allq_Q(2*b,h,H,g,1,LSflag,ndates,Euler);
sol = sol.'.*exp(-alpha*x);
price = price_LB(sol,x,lx,lk,param);
CPU_Time = toc;
Ris(1,1:2)=[price, CPU_Time];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('- Quadrature with 4th order method (as Simpson)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
sol = solve_allq_Q(2*b,h,H,g,4,LSflag,ndates,Euler);
sol = sol.'.*exp(-alpha*x);
price = price_LB(sol,x,lx,lk,param);
CPU_Time = toc;
Ris(1,3:4)=[price, CPU_Time];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('- ZT + FHT with sinc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
sol = solve_allq_FHT_DB(H,G,l,u,b,ndates,Euler,@(x)ifht(x,'sinc',0));
sol=real(sol(end/2+1:end)).*exp(-alpha*x);   
price = price_LB(sol,x,lx,lk,param);
CPU_Time = toc;
S=exp(lx-x); S=fliplr(S);
cdfmin=1-sol; cdfmin=fliplr(cdfmin);
cdfmin=max(cdfmin,0);
[m,index]=min(cdfmin); cdfmin(1:index-1)=m*(S(1:index-1)-S(1))/(S(index)-S(1));
SS=linspace(S(1),K,ngrid);
price2=exp(-rf*T)*trapz(SS,interp1(S,cdfmin,SS,'spline'));
Ris(1,5:6)=[price, CPU_Time];
% figure
% plot(S,cdfmin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Spitzer with sinc')
% BASATO SU Green, Fusai, Abrahams (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=S_0;
alpha = 0; % shift parameter = damping factor

b = up; % upper bound of the support
trunc_L = S_0*exp(-b); % lower truncation
trunc_U = S_0*exp(b); % upper truncation

[x,h,xi,H] = kernel(ngrid,-b,b,param,alpha,1,1);
[S,g,G] = payoff(x,xi,alpha,K,trunc_L,trunc_U,C,-1);
%[S,g,G] = payoff_old(x,xi,alpha,K,trunc_L,trunc_U,C,-1);

tic
sol = solve_allq_FHT_S_LB(H,G,b,ndates,Euler,@(x)ifht(x,'sinc',ngrid),-1,x,C,S_0);
price = sol;
CPU_Time = toc;
Ris(1,7:8)=[price, CPU_Time];

fprintf('%.6f - %.2f -- %.6f - %.2f -- %.10f - %.2f -- %.10f - %.2f\n',Ris);
end