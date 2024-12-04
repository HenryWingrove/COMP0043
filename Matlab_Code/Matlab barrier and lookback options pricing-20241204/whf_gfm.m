% Solve the Wiener-Hopf (d == -b or u == b) or Fredholm (-b < d < u < b)
% equation f(x) = \int_d^u f(y)h(x-y) dy + g(x) by fast Hilbert transforms.
% The kernel h and the forcing function g are given on the support [-b,b]
% with b <= d < u <= b, and f is unknown. The integral of the kernel h,
% which is equal to H(0), must be < 1, otherwise log(1-H) diverges. Here
% several kernels and forcing functions can be given as rows of the input
% matrices H and G; the solutions are the rows of the returned matrix F.
% Annotated references:
% 1. R. J. Henery, J. Inst. Maths. Applics. 13, 89-96 (1974):
%    the steps (i)-(iii) on page 93 solve the Wiener-Hopf equation,
%    though he does not recognise the Hilbert transform.
% 2. J. A. C. Weideman, Math. Comp. 64, 745-762 (1995):
%    four ways to compute numerically the Hilbert transform.
% 3. G. Vergara-Caffarelli, P. Loreti, Trasformata di Hilbert,
%    Università  di Roma La Sapienza, 1999: Plemelj-Sokhotski formulas.
% 4. R. J. Henery, J. Math. Phys. ? (1981): how to compute the Hilbert
%    transform of a diverging argument, which happens here if H(0) >= 1.
% 5. R. Green, G. Fusai, I. D. Abrahams, Math. Finance 20, 259-288 (2010):
%    Section 2.4, for the Fredholm case.
% Authors: Guido Germano, Daniele Marazzina, Gianluca Fusai, 2009-2012.
function [F0,miter] = whf_gfm(q,H,G,d,u,b,ifht)

% Build input matrices in Fourier space
%H = fftshift(ifft(ifftshift(h,2),[],2),2)*2*b; % characteristic function
%G = fftshift(ifft(ifftshift(g,2),[],2),2)*2*b;
m = length(q);
n = length(H);
H = repmat(q,1,n).*repmat(H,m,1);
G = repmat(G,m,1);
xi = pi/b*(-n/2:n/2-1); % grid in Fourier space
miter = 0;

% Factorise L = 1/(1-H) with respect to zero
lL = -log(1-H); % 1-h is the functional derivative of the eq. w.r.t. f
iHlL = ifht(lL); % imaginary unit times the fast Hilbert transform of L
lLp = (lL+iHlL)/2; % Plemelj-Sokhotski
lLm = (lL-iHlL)/2; % Plemelj-Sokhotski
Lp = exp(lLp);
Lm = exp(lLm);

if -b < d && u == b % Wiener-Hopf equation on (d,+\infty)

    % Decompose Cd = G.*Ldm with respect to the down barrier
    eidxi = repmat(exp(1i*d*xi),m,1);
    Cd = G.*Lm;
    ifHCd = eidxi.*ifht(Cd./eidxi); % 1i*FHT of Cd with respect to d
    Cdp = (Cd+ifHCd)/2; % Plemelj-Sokhotski
    
    % Compute the solution
    F0 = Cdp.*Lp;
    miter = miter + 1;
    
elseif -b == d && u < b % Wiener-Hopf equation on (-\infty,u)
    
    % Decompose Cu = G.*Lup with respect to the up barrier
    eiuxi = repmat(exp(1i*u*xi),m,1);
    Cu = G.*Lp;
    ifHCu = eiuxi.*ifht(Cu./eiuxi); % 1i*FHT of Cu with respect to u
    Cum = (Cu-ifHCu)/2; % Plemelj-Sokhotski
    
    % Compute the solution
    F0 = Cum.*Lm;
    miter = miter + 1;
    
else % Fredholm equation on (d,u)
    
    maxiter = 5;
    tol = 1e-12;
    index = n/4:3*n/4;
    
    eidxi = exp(1i*d*xi);
    eiuxi = exp(1i*u*xi);
    F0 = zeros(m,n);
    for i = 1:m
        Fu = zeros(1,n); % initial guess
        iter = 0;
        while 1
            
            % Decompose Cd = (G-Fu).*Lm with respect to the down barrier
            Cd = (G(i,:)-Fu).*Lm(i,:);
            ifHCd = eidxi.*ifht(Cd./eidxi); % if d == -b, ifHCd = Cd
            Cdp = (Cd+ifHCd)/2; % Plemelj-Sokhotski
            
            % Compute the solution and its diff. w.r.t. the previous iteration
            F0_new = Cdp.*Lp(i,:);
            err = norm(F0_new(index)-F0(i,index),'inf');
            F0(i,:) = F0_new;
            iter = iter + 1;
            if err < tol || iter == maxiter
                break;
            end
            
            % Prepare the next iteration
            Cdm = (Cd-ifHCd)/2; % Plemelj-Sokhotski
            Fd = Cdm./Lm(i,:);
            
            
            % Decompose Cu = (G-Fd).*Lp with respect to the up barrier
            Cu = (G(i,:)-Fd).*Lp(i,:);
            ifHCu = eiuxi.*ifht(Cu./eiuxi); % if u == b, ifHCu = -Cu
            Cum = (Cu-ifHCu)/2; % Plemelj-Sokhotski
            
            % Compute the solution and its diff. w.r.t. the previous iteration
            F0_new = Cum.*Lm(i,:);
            err = norm(F0_new(index)-F0(i,index),'inf');
            F0(i,:) = F0_new;
            iter = iter + 1;
            if err < tol || iter == maxiter
                break;
            end
            
            % Prepare the next iteration
            Cup = (Cu+ifHCu)/2; % Plemelj-Sokhotski
            Fu = Cup./Lp(i,:);

        end
        miter = miter + iter;

        %if i == 1
        %    figure(1), clf
        %    x = b/(2*n)*(-n/2:n/2-1);
        %    f0 = real(fftshift(fft(ifftshift(F0_new))))/(2*b);
        %    fd = real(fftshift(fft(ifftshift(Fd))))/(2*b);
        %    fu = real(fftshift(fft(ifftshift(Fu))))/(2*b);
        %    plot(x,f0,'r',x,fd,'g',x,fu,'b')
        %    legend('f0','fd','fu')
        %end

    end
    miter = miter/m;

end

%f0 = real(fftshift(fft(ifftshift(F0,2),[],2),2))/(2*b);
