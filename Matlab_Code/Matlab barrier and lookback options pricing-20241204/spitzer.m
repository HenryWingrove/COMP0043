% See R. Green, G. Fusai, I. D. Abrahams, Math. Finance 20, 259-288 (2010),
% Section 3.1.
% Authors: Guido Germano, Daniele Marazzina, Gianluca Fusai, 2012-2014
function F0 = spitzer(q,K,G,d,b,ifht)

% Build input matrices in Fourier space
%H = fftshift(ifft(ifftshift(h,2),[],2),2)*2*b; % characteristic function
m = length(q);
n = length(K);
K = repmat(K,m,1);
L = 1-repmat(q,1,n).*K;
%G = repmat(G,m,1);
xi = pi/b*(-n/2:n/2-1); % grid in Fourier space

% Factorise L = (1-qK) with respect to zero
lL = log(L); % GFA (2.11)
iHlL = ifht(lL); % imaginary unit times the fast Hilbert transform of L
lLm = (lL-iHlL)/2; % Plemelj-Sokhotski
Lm = exp(lLm);
lLp = (lL+iHlL)/2; % Plemelj-Sokhotski
Lp = exp(lLp);

% Decompose P with respect to the down barrier
eidxi = repmat(exp(1i*d*xi),m,1);
P = K./(eidxi.*Lm); % GFA (2.12)
ifHP = ifht(P); % 1i*FHT of P with respect to 0
Pp = (P+ifHP)/2; % Plemelj-Sokhotski

% Compute the solution
%F0 = K.*conj(G).*(K-eidxi.*Pm.*Lm)./L; % GFA (3.6.c)
%F0 = K.*conj(G).*eidxi.*Pp./Lp; % GFA (3.6.c)
F0 = eidxi.*Pp./Lp; % GFA (3.6.c) %% SE MOLTIPLICO PER K*conj(G) FUORI!!!