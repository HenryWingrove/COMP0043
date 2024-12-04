% See R. Green, G. Fusai, I. D. Abrahams, Math. Finance 20, 259-288 (2010),
% Section 3.1.
% Authors: Guido Germano, Daniele Marazzina, Gianluca Fusai, June 2012.
function F0 = spitzer_LB(q,H,ifht,call)

% Build input matrices in Fourier space
%H = fftshift(ifft(ifftshift(h,2),[],2),2)*2*b; % characteristic function
m = length(q);
n = length(H);
H = repmat(q,1,n).*repmat(H,m,1);

% Factorise L = 1/(1-H) with respect to zero
lL = -log(1-H); % GFA (2.11)
iHlL = ifht(lL); % imaginary unit times the fast Hilbert transform of L
lLm = (lL-iHlL)/2; % Plemelj-Sokhotski
lLp = (lL+iHlL)/2; % Plemelj-Sokhotski
Lm = exp(lLm);
Lp = exp(lLp);

if call == 1 % call su massimo
    F0 = repmat(Lm(:,end/2+1),1,n).*Lp;
else % put su minimo
    F0 = repmat(Lp(:,end/2+1),1,n).*Lm;
end