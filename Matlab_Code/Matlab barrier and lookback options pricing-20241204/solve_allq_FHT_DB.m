function [sol,miter] = solve_allq_FHT_DB(H,Sol,d,u,b,ndates,Euler,ifht,acr)

if nargin==8
    acr=0; %attiva l'accrocchio per avere comportamento esponenziale
end
miter=0;
if ndates > 1

    gamma = 6; % accuracy = 10^(-2*gamma)
    n = 12; m = 20; % parameters for Euler acceleration/summation

    if acr == 0
        if ndates == 2
            ndates = 1;
        else
            Sol = Sol.*H;
            ndates = ndates - 2;
        end
    end
    
    rho = 10^(-gamma/ndates);

    if ndates <= n+m || ~Euler % Standard method
        % G. Fusai, I. D. Abrahams, C. Sgarra
        % Finance Stochast. 10, 1-26 (2006), Eq. (26)
        q = rho*[1;exp(-1i*(1:ndates-1).'*pi/ndates);-1];
        c = (-1).^(0:ndates).'.*[0.5;ones(ndates-1,1);0.5];
        [temp,miter] = whf_gfm(q,H,Sol,d,u,b,ifht);
%         C=-2;
%         if abs(C)>0
%             for iii=1:length(temp)
%                 temp(:,iii)=temp(:,iii).*(q.^-C);
%             end
%         end
        Sol = sum(repmat(c,size(H)).*temp,1);
    else % Euler acceleration/summation
        % Approximate with the binomial average or Euler transform
        % of the partial sums of order n+1 to n+m
        q = rho*[1;exp(-1i*(1:n+m).'*pi/ndates)];
        c = (-1).^(0:n+m).'.*[0.5;ones(n+m,1)];
        [temp,miter] = whf_gfm(q,H,Sol,d,u,b,ifht);
        ps = cumsum(repmat(c,size(H)).*temp,1);
        Sol = sum(repmat(binomial(m).',size(H)).*ps(n:n+m,:),1)/2^m;
        %figure, clf, hold on
        %plot(real(ps(:,end/2+10)),'or')
        %plot(real(repmat(Sol(end/2+10),size(ps,1))),'r')
        %plot(imag(ps(:,end/2+10)),'sg')
        %plot(imag(repmat(Sol(end/2+10),size(ps,1))),'g')
    end
    Sol = Sol/(ndates*rho^ndates);

end
if acr==0
    Sol = Sol.*H;
end
sol = real(fftshift(fft(ifftshift(Sol))))/(2*b);