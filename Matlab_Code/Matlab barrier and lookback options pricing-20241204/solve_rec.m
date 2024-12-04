function sol = solve_rec(H,g,Qflag,ndates)

% Parameters
N = length(g)/2;

% Quadrature weights
w = weights(Qflag,N);

sol = transpose(g(N+1:end));
H = ifftshift(transpose(H));
for i = 1:ndates
    sol = product(sol,H,w);
end
    
function z = product(r,H,w)

N = length(r);
U = ifft([r.*w;zeros(N,1)]);
v = fft(H.*U);
z = real(v(1:N));
