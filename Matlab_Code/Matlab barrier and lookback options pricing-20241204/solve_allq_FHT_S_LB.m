function sol = solve_allq_FHT_S_LB(H,G,b,ndates,Euler,ifht,call,x,C,S_0)

gamma = 6; % accuracy = 10^(-2*gamma)
n = 12; m = 20; % parameters for Euler acceleration/summation
rho = 10^(-gamma/ndates);

ndates = ndates-1;

if ndates <= n+m || ~Euler % Standard method
    % G. Fusai, I. D. Abrahams, C. Sgarra
    % Finance Stochast. 10, 1-26 (2006), Eq. (26)
    q = rho*[1;exp(-1i*(1:ndates-1).'*pi/ndates);-1];
    c = (-1).^(0:ndates).'.*[0.5;ones(ndates-1,1);0.5];
    temp = spitzer_LB(q,H,ifht,call);
    Sol = sum(repmat(c,size(H)).*temp,1);
else % Euler acceleration/summation
    % Approximate with the binomial average or Euler transform
    % of the partial sums of order n+1 to n+m
    q = rho*[1;exp(-1i*(1:n+m).'*pi/ndates)];
    c = (-1).^(0:n+m).'.*[0.5;ones(n+m,1)];
    temp = spitzer_LB(q,H,ifht,call);
    ps = cumsum(repmat(c,size(H)).*temp,1);
    Sol = sum(repmat(binomial(m).',size(H)).*ps(n:n+m,:),1)/2^m;
    %figure, clf, hold on
    %plot(real(ps(:,end/2+10)),'or')
    %plot(real(repmat(Sol(end/2+10),size(ps,1))),'r')
    %plot(imag(ps(:,end/2+10)),'sg')
    %plot(imag(repmat(Sol(end/2+10),size(ps,1))),'g')
end
Sol = Sol/(ndates*rho^ndates);
Sol = Sol.*H;
sol = real(fftshift(fft(ifftshift((Sol).*conj(G)))))/(2*b);
if C==S_0
    sol = sol(end/2+1); % x=0
else
    sol=interp1(x,sol,log(S_0/C),'spline');
end