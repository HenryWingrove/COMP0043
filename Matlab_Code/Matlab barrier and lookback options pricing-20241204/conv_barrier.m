% Fourier Space Time-stepping method for pricing barrier options.
% See Ken Jackson, Sebastian Jaimungal, Vladimir Surkov,
% "Fourier space time-stepping for option pricing with Levy models",
% http://papers.ssrn.com/abstract=1020209 (2007).
% Gianluca Fusai and Guido Germano, September 2010.

function sol = conv_barrier(H,sol,ind_fun,ndates)

%fftw('planner','measure');
H = ifftshift(H);
for j = 1:ndates
    sol = real(fft(ifft(sol).*H));
    sol(ind_fun) = 0;
end
