% Solve the Wiener-Hopf equation f(t) = \int_0^\infty h(t-s) f(s) ds + g(t)
% by composite Newton-Cotes quadrature. The kernel h and the forcing function g
% are given on a support with length T symmetric around 0, and f is unknown.
function [f,iter] = quadrature(T,h,H,g,Qflag,LSflag)

% Parameters
DT = T/length(h);
N = length(h)/2;

% Quadrature weights
w = weights(Qflag,N);

% Forcing function
g = transpose(g(N+1:end));

if LSflag == 0 % Gaussian elimination, O(N^3)
    hmat = toeplitz(h(N+1:end),h(N+1:-1:2)); % column, row
    amat = eye(N) - hmat*DT.*repmat(transpose(w),N,1);
    f = amat\g;
    iter = 1;
else % Generalized minimum residual method, O(N*log(N))
    tol = 1e-12;
    maxit = 100;
    H = ifftshift(transpose(H));
    if LSflag == 1 % No preconditioner
        [f,dummy,dummy,iter] = gmres(@product,g,[],tol,maxit,[],[],g,H,[],w);
    else % Lin, Ng and Chan preconditioner
        P = H./(H-1); % -H contains the eigenvalues of the circulant matrix
        [f,dummy,dummy,iter] = gmres(@product,g,[],tol,maxit,@precond,[],g,H,P,w);
    end
    f = real(f);
    iter = iter(2);
end

% Compute the product p = (eye(N)-hmat*DT)*r, where hmat is a Toeplitz
% matrix built from h, through a circulant matrix specified by H and fast FTs,
% exploiting that hmat*DT*r = v(N+1:end), where V = H.*U and u = [zeros(N,1);r].
% The factor DT = T/(2*N) simplifies away because a factor T is already
% included in H and a factor 2*N comes from v = 2*N*fft(H.*U).
% The factor 2*N is not necessary if fft and ifft are swapped.
function z = product(r,H,dummy,w)

U = ifft([r.*w;zeros(size(r))]);
v = fft(H.*U);
z = r - v(1:end/2);

% Preconditioner: FR Lin, MK Ng, RH Chan, SIAM J. Numer. Anal. 34, 1418 (1997)
function z = precond(r,dummy,P,w)

U = ifft([r.*w;zeros(size(r))]);
v = fft(P.*U);
z = r - v(1:end/2);
