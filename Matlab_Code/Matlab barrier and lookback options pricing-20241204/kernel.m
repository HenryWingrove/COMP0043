function [x,h,w,H] = kernel(ngrid,xmin,xmax,parameters,alpha,disc,flag)
% disc=0 --> no discount factor in the density
% disc=1 --> discount factor in the density
% flag=0 --> funzione caratteristica per problema backward 
% flag=1 --> funzione caratteristica per problema forward

if nargin==5
    disc=1; %discount factor
    flag=0; %backward characteristic function
elseif nargin==6
    flag=0; %backward characteristic function
end
    
N = ngrid/2;
dx = (xmax-xmin)/ngrid;
x = dx*(-N:N-1);
dw = 2*pi/(xmax-xmin);
w = dw*(-N:N-1);
if nargin < 5 % shift parameter, esp. for Feng-Linetsky and convolution
    alpha = 0;
end
H = charfunction(w+1i*alpha,parameters,flag); % characteristic function
if disc==1
    H = H*exp(-parameters.rf*parameters.dt); % discount
end
h = real(fftshift(fft(ifftshift(H))))/(xmax-xmin); % discounted kernel


