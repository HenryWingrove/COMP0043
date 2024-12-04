function sol = feng_linetsky(H,Sol,d,u,b,ndates,ifht)
% L. Feng, V. Linetsky, Math. Finance 18, 337-384 (2008)

n = length(H);
xi = pi/b*(-n/2:n/2-1); % grid in Fourier space

if -b < d && u == b % down-and-out option
    eidxi = exp(1i*d*xi);
    for i = ndates:-1:2 % backward induction in Fourier space
        Sol = Sol.*H;
        Sol = 0.5*(Sol+eidxi.*ifht(Sol./eidxi)); % Eq. (5.11)
    end
elseif -b == d && u < b % up-and-out option
    eiuxi = exp(1i*u*xi);
    for i = ndates:-1:2 % backward induction in Fourier space
        Sol = Sol.*H;
        Sol = 0.5*(Sol-eiuxi.*ifht(Sol./eiuxi)); % Eq. (5.12)
    end
else % double-barrier option
    eidxi = exp(1i*d*xi);
    eiuxi = exp(1i*u*xi);
   %sp = (u+d)/2;
   %sm = (u-d)/2;
   %esvec = exp(1i*sp*xi).*sin(sm*xi)./(pi*xi);
   %esvec(end/2+1) = sm/pi;
   %esvec = real(fft(ifftshift(esvec)));
    for i = ndates:-1:2 % backward induction in Fourier space
        Sol = Sol.*H;
        Sol = 0.5*(eidxi.*ifht(Sol./eidxi)-eiuxi.*ifht(Sol./eiuxi)); % Eq. (5.10)
       %Sol = fft(esvec.*ifft(Sol)); % Eq. (5.13)
    end
end

% Eqs. (5.14) and (6.38): final step, from Fourier space to normal space
Sol = Sol.*H;
sol = real(fftshift(fft(ifftshift(Sol))))/(2*b);