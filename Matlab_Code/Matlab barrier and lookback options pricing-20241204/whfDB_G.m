% See R. Green, G. Fusai, I. D. Abrahams, Math. Finance 20, 259-288 (2010),
% Section 2.4.
function [F0,miter] = whfDB_G(q,K,G,d,u,b,ifht)

% Build input matrices in Fourier space
%H = fftshift(ifft(ifftshift(h,2),[],2),2)*2*b; % characteristic function
m = length(q);
n = length(K);
L = 1-repmat(q,1,n).*repmat(K,m,1);
xi = pi/b*(-n/2:n/2-1); % grid in Fourier space

% Factorise L = 1/(1-q.*K) with respect to zero
lL = -log(L); % GFA (2.11)
iHlL = ifht(lL); % imaginary unit times the fast Hilbert transform of L
lLm = (lL-iHlL)/2; % Plemelj-Sokhotski
Lm = exp(lLm);
lLp = (lL+iHlL)/2; % Plemelj-Sokhotski
Lp = exp(lLp);

eidxi = exp(1i*d*xi);
eiuxi = exp(1i*u*xi);
KG = K.*G;

miter = 0;
maxiter = 5;
tol = 1e-10;
F0 = zeros(m,n);
for i = 1:m
        Fu = zeros(1,n); % initial guess
        iter = 0;
        while 1
            
            % Decompose Cd = (KG-Fu).*Lm with respect to the down barrier
            Cd = (KG-Fu).*Lm(i,:);
            ifHCd = eidxi.*ifht(Cd./eidxi); % if d == -b, ifHCd = Cd
            Cdp = (Cd+ifHCd)/2; % Plemelj-Sokhotski
            
            % Compute the solution and its diff. w.r.t. the previous iteration
            F0_old = F0(i,:);
            F0(i,:) = Cdp.*Lp(i,:);
            err = norm(F0(i,:)-F0_old,'inf');
            iter = iter + 1;
            if err < tol || iter == maxiter
                break;
            end
            
            % Prepare the next iteration
            Cdm = (Cd-ifHCd)/2; % Plemelj-Sokhotski
            Fd = Cdm./Lm(i,:);
                        
            % Decompose Cu = (KG-Fd).*Lp with respect to the up barrier
            Cu = (KG-Fd).*Lp(i,:);
            ifHCu = eiuxi.*ifht(Cu./eiuxi); % if u == b, ifHCu = -Cu
            Cum = (Cu-ifHCu)/2; % Plemelj-Sokhotski
            
            % Compute the solution and its diff. w.r.t. the previous iteration
            F0_old = F0(i,:);
            F0(i,:) = Cum.*Lm(i,:);
            err = norm(F0(i,:)-F0_old,'inf');
            iter = iter + 1;
            if err < tol || iter == maxiter
                break;
            end
            
            % Prepare the next iteration
            Cup = (Cu+ifHCu)/2; % Plemelj-Sokhotski
            Fu = Cup./Lp(i,:);

        end
        miter = miter + iter;
end
miter = miter/m;