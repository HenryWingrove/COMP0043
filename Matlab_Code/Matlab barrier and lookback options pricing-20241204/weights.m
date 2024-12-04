function w = weights(Qflag,N)
% Quadrature weights

if Qflag == 0 % rectangles
    w = ones(N,1);
elseif Qflag == 1 % trapezia, O(1/N^2), requires N > 1
    w = [0.5;ones(N-2,1);0.5];
elseif Qflag == 2 % Simpson, O(1/N^4), requires odd N > 2
    w = ones(N,1);
    w(2:2:end) = 4;
    w(3:2:end-1) = 2;
    w = w/3;
elseif Qflag == 3 % 3rd order method, Eq. (4.1.12) Numerical Recipes 3rd ed.
    w = [5/12;13/12;ones(N-4,1);13/12;5/12];
elseif Qflag == 4 % other 4th order method, Eq. (4.1.14) Numerical Recipes 3rd ed.
    w = [3/8;7/6;23/24;ones(N-6,1);23/24;7/6;3/8];
end
