function [sol,iter_vec,iter_ave,iter_max] = solve_allq_Q(T,h,H,g,Qflag,LSflag,ndates,Euler)

% G. Fusai, I. D. Abrahams, C. Sgarra
% Finance Stochast. 10, 1-26 (2006), Eq. (26)

gamma = 6;
rho = 10^(-gamma/ndates);
m = 20; n = 12; % Parameters for Euler acceleration/summation

iter_ave = 0;
iter_max = 0;
iter_vec = 0*(1:min(ndates+1,m+n+1));
if ndates <= m+n || ~Euler % Standard method
    q = rho*[1 exp(-1i*(1:ndates-1)*pi/ndates) -1];
    coeff = [1 2*ones(1,ndates-1) 1];
    sol = 0;
    for j = 1:ndates+1
        [temp,iter] = quadrature(T,q(j)*h,q(j)*H,g,Qflag,LSflag);
        temp = real(temp);
        sol = sol + (-1)^(j-1)*coeff(j)*temp;
        iter_vec(j) = iter;
        iter_ave = iter_ave + iter;
        iter_max = max(iter_max,iter);
    end
    sol = sol/(2*ndates*rho^ndates);
    iter_ave = iter_ave/(ndates+1);
else % Euler acceleration/summation
    q = rho*[1 exp(-1i*(1:m+n)*pi/ndates)];
    coeff = [0.5 ones(1,n)];
    SN = 0;
    for j = 1:n+1
        [temp,iter] = quadrature(T,q(j)*h,q(j)*H,g,Qflag,LSflag);
        temp = real(temp);
        SN = SN + (-1)^(j-1)*coeff(j)*temp;
        iter_vec(j) = iter;
        iter_ave = iter_ave + iter;
        iter_max = max(iter_max,iter);
    end
    sol = binomial(m,0)*SN;
    for j = 1:m
        [temp,iter] = quadrature(T,q(j+n+1)*h,q(j+n+1)*H,g,Qflag,LSflag);
        temp = real(temp);
        SN = SN + (-1)^j*temp;
        sol = sol + binomial(m,j)*SN;
        iter_vec(j+n+1) = iter;
        iter_ave = iter_ave + iter;
        iter_max = max(iter_max,iter);
    end
    sol = sol/(2^m*rho^ndates*ndates);
    iter_ave = iter_ave/(m+n+1);
end
