% See R. Green, G. Fusai, I. D. Abrahams, Math. Finance 20, 259-288 (2010),
% Section 2.4.
function [F0,miter] = whfDB_D(q,K,G,d,u,b,ifht)

% Build input matrices in Fourier space
%H = fftshift(ifft(ifftshift(h,2),[],2),2)*2*b; % characteristic function
m = length(q);
n = length(K);
L = 1-repmat(q,1,n).*repmat(K,m,1);
xi = pi/b*(-n/2:n/2-1); % grid in Fourier space

% Factorise L = 1-qK with respect to zero
lL = log(L); % GFA (2.11)
iHlL = ifht(lL); % imaginary unit times the fast Hilbert transform of L
lLm = (lL-iHlL)/2; % Plemelj-Sokhotski
Lm = exp(lLm);
lLp = (lL+iHlL)/2; % Plemelj-Sokhotski
Lp = exp(lLp);

eidxi = exp(1i*d*xi);
eiuxi = exp(1i*u*xi);
emidxi = exp(-1i*d*xi);
emiuxi = exp(-1i*u*xi);
eiudxi = exp((u-d)*1i*xi);
eiduxi = exp((d-u)*1i*xi);

miter = 0;
maxiter = 5;
tol = 1e-10;
F0 = zeros(m,n);
for i = 1:m
    Jp = zeros(1,n); % initial guess
    Jm = zeros(1,n); % initial guess
    iter = 0;
    while 1        
        P = (emidxi.*K.*G-eiudxi.*Jp)./Lm(i,:);
        ifHP = ifht(P); % 1i*FHT of P with respect to 0
        Pm = (P-ifHP)/2; % Plemelj-Sokhotski

        Q = (emiuxi.*K.*G-eiduxi.*Jm)./Lp(i,:);
        ifHQ = ifht(Q); % 1i*FHT of Q with respect to 0
        Qp = (Q+ifHQ)/2; % Plemelj-Sokhotski

        Jm = Pm.*Lm(i,:); % GFA (2.55)
        Jp = Qp.*Lp(i,:); % GFA (2.56)

        % Compute the solution and its diff. w.r.t. the previous iteration
        F0_old = F0(i,:);
        F0(i,:) = K.*(K.*G-eidxi.*Jm-eiuxi.*Jp)./L(i,:); % GFA (2.60)
        %F0(i,:) = K.*G.*(K./L(i,:)-eidxi.*Pm./Lp(i,:)-eiuxi.*Qp./Lm(i,:)); % GFA (2.60)           
        err = norm(F0(i,:)-F0_old,'inf');
        iter = iter + 1;
        if err < tol || iter == maxiter
            break;
        end
    end
    miter = miter + iter;
end
miter = miter/m;